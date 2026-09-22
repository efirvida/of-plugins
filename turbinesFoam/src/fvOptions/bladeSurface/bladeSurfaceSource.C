/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam, which is based on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "bladeSurfaceSource.H"
#include "actuatorLineElement.H"
#include "volFields.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::bladeSurfaceSource::createOutputFiles()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path()/"../postProcessing/bladeSurface";
    }
    else
    {
        dir = mesh_.time().path()/"postProcessing/bladeSurface";
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    stationFile_ = new OFstream(dir/name_ + ".csv");

    // The per-node CSV is the source of the design's surface-moment
    // reconstruction oracle (rtol 1e-6, design section 8.1); the default six
    // significant digits cannot represent it. Per-stream precision, so the
    // element and turbine CSVs keep the delivered precision (and the ALM /
    // no-mesh ASM outputs stay byte-identical).
    stationFile_->precision(12);

    *stationFile_
        << "time,station,root_dist,area,force_x,force_y,force_z,"
        << "c_ref_n,c_ref_t,f_ref_n,f_ref_t" << endl;

    if (writeNodePerf_)
    {
        nodeFile_ = new OFstream(dir/name_ + "_nodes.csv");
        nodeFile_->precision(12);

        *nodeFile_
            << "time,node,x,y,z,nx,ny,nz,fx,fy,fz,area,station,"
            << "chord_fraction" << endl;
    }
}


void Foam::fv::bladeSurfaceSource::writeOutput()
{
    if (writePerf_ and Pstream::master())
    {
        writeStationCsv();

        if (writeNodePerf_)
        {
            writeNodeCsv();
        }
    }
}


void Foam::fv::bladeSurfaceSource::writeStationCsv()
{
    // Patch force on the blade, per unit density: partition of unity makes
    // the patch sum equal the element force
    List<vector> patchForce(nElements_, vector::zero);

    forAll(sampler_.nodeForces_, i)
    {
        patchForce[sampler_.patch_[i]] += sampler_.nodeForces_[i];
    }

    const scalar rhoRef = sampler_.rhoRef();
    const scalar time = mesh_.time().value();

    forAll(elements_, e)
    {
        const vector fSI = rhoRef*patchForce[e];

        // One element per station, so the element's public reference
        // quantities are the station values (the same definitions the element
        // CSV writes, D9); root_dist uses the element convention that
        // comparePhaseVI.py maps to r/R
        *stationFile_ << time
            << "," << sampler_.elementStation()[e]
            << "," << elements_[e].rootDistance()
            << "," << sampler_.patchArea()[e]
            << "," << fSI.x() << "," << fSI.y() << "," << fSI.z()
            << "," << elements_[e].normalRefCoefficient()
            << "," << elements_[e].tangentialRefCoefficient()
            << "," << elements_[e].normalRefForce()
            << "," << elements_[e].tangentialRefForce()
            << endl;
    }
}


void Foam::fv::bladeSurfaceSource::writeNodeCsv()
{
    // Body-frame positions/normals and SI on-blade forces (the contract)
    const List<point>& X = sampler_.positions();
    const List<vector>& n = sampler_.normals();
    const List<vector>& F = sampler_.forces();
    const List<scalar>& A = sampler_.areas();
    const List<scalar>& s = sampler_.station();
    const List<scalar>& c = sampler_.chordFraction();
    const scalar time = mesh_.time().value();

    forAll(F, i)
    {
        *nodeFile_ << time << "," << i
            << "," << X[i].x() << "," << X[i].y() << "," << X[i].z()
            << "," << n[i].x() << "," << n[i].y() << "," << n[i].z()
            << "," << F[i].x() << "," << F[i].y() << "," << F[i].z()
            << "," << A[i]
            << "," << s[i]
            << "," << c[i] << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::bladeSurfaceSource::bladeSurfaceSource
(
    const word& ownerName,
    const dictionary& dict,
    const fvMesh& mesh,
    PtrList<actuatorLineElement>& elements
)
:
    sampler_(dict, mesh, elements),
    mesh_(mesh),
    elements_(elements),
    nElements_(elements.size()),
    name_(ownerName + ".surface"),
    writePerf_(dict.lookupOrDefault("writePerf", true)),
    writeNodePerf_(dict.lookupOrDefault("writeNodePerf", false)),
    stationFile_(nullptr),
    nodeFile_(nullptr)
{
    if (writePerf_)
    {
        createOutputFiles();
    }

    Info<< "Blade surface source '" << name_ << "': "
        << sampler_.nNodes() << " nodes, " << nElements_
        << " element patches" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::bladeSurfaceSource::~bladeSurfaceSource()
{
    delete stationFile_;
    delete nodeFile_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::bladeSurfaceSource::distribute
(
    volVectorField& bladeForceField,
    vector& bladeForce,
    const volScalarField* rhoPtr
)
{
    // Compressible overload: the element's public force() carries the density
    // at the element position (multiplyForceRho), so the element-local density
    // is needed to recover the per-unit-density share
    List<scalar> rhoElement(nElements_, 1.0);

    if (rhoPtr)
    {
        forAll(elements_, e)
        {
            scalar rhoE = VGREAT;

            const label cellI = mesh_.findCell(elements_[e].position());

            if (cellI >= 0)
            {
                rhoE = (*rhoPtr)[cellI];
            }

            reduce(rhoE, minOp<scalar>());

            if (!(rhoE < VGREAT))
            {
                FatalErrorInFunction
                    << "Blade surface element " << e << " at "
                    << elements_[e].position() << " not found in mesh" << nl
                    << exit(FatalError);
            }

            rhoElement[e] = rhoE;
        }
    }

    // Per-node share of the element force, per unit density (D4.3), and the
    // force actually distributed (node-local density in the compressible
    // overload)
    List<vector> distributed(sampler_.nNodes(), vector::zero);

    forAll(sampler_.nodeForces_, i)
    {
        const label e = sampler_.patch_[i];

        vector F = sampler_.patchAreaShare_[i]*elements_[e].force();

        if (rhoPtr)
        {
            F /= max(rhoElement[e], VSMALL);
        }

        sampler_.nodeForces_[i] = F;
        sampler_.forces_[i] = sampler_.rhoRef()*(sampler_.bodyToGlobal().T() & F);

        if (rhoPtr)
        {
            scalar rhoNode = VGREAT;

            const label cellI =
                mesh_.findCell(sampler_.positionsGlobal()[i]);

            if (cellI >= 0)
            {
                rhoNode = (*rhoPtr)[cellI];
            }

            reduce(rhoNode, minOp<scalar>());

            if (!(rhoNode < VGREAT))
            {
                FatalErrorInFunction
                    << "Blade surface sample at " << sampler_.positionsGlobal()[i]
                    << " not found in mesh" << nl << exit(FatalError);
            }

            distributed[i] = F*rhoNode;
        }
        else
        {
            distributed[i] = F;
        }
    }

    // Eq. 18 over the bounded candidate lists
    sampler_.distributeNodeForces(bladeForceField, distributed);

    // Distributed total (global). It replaces the element-loop total so the
    // reported blade force matches the applied load; in the incompressible
    // reference partition of unity makes it equal the summed element forces
    vector total = vector::zero;

    forAll(distributed, i)
    {
        total += distributed[i];
    }

    bladeForce = returnReduce(total, sumOp<vector>());

    writeOutput();
}


Foam::vector Foam::fv::bladeSurfaceSource::moment(const point& p) const
{
    vector m = vector::zero;
    const List<point>& X = sampler_.positionsGlobal();

    forAll(X, i)
    {
        m += (X[i] - p) ^ sampler_.nodeForces_[i];
    }

    return m;
}


// ************************************************************************* //
