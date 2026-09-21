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

#include "nacelleSurfaceSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "volFields.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(nacelleSurfaceSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        nacelleSurfaceSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::nacelleSurfaceSource::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path()/"../postProcessing/nacelle";
    }
    else
    {
        dir = mesh_.time().path()/"postProcessing/nacelle";
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_ = new OFstream(dir/name_ + ".csv");

    *outputFile_ << "time,fx,fy,fz,f_n_mag,f_tau_mag,cd" << endl;

    if (writeNodePerf_)
    {
        nodeFile_ = new OFstream(dir/name_ + "_nodes.csv");

        *nodeFile_
            << "time,node,x,y,z,nx,ny,nz,fx,fy,fz,area" << endl;
    }
}


void Foam::fv::nacelleSurfaceSource::calcForceField
(
    const volVectorField& U,
    const volScalarField* rhoPtr
)
{
    const scalar dt = mesh_.time().deltaT().value();

    // Zero the force field and the totals
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    force_ = vector::zero;
    forceNormal_ = vector::zero;
    forceTangential_ = vector::zero;

    const List<point>& X = sampler_.positionsGlobal();
    const List<scalar>& A = sampler_.areas();

    List<scalar> nodeH(X.size(), 0.0);

    forAll(X, i)
    {
        // Local cell size h = cbrt(V[cell]); fatal on an unreachable sample
        const scalar h = sampler_.cellSize(X[i]);
        nodeH[i] = h;

        // Eq. 7: kernel-smoothed velocity at the node
        const vector uTilde = sampler_.interpolateVelocity(X[i], U, h);

        // Eq. 19: normal traction (u^d = 0)
        const vector fN = sampler_.normalForce(i, uTilde, h, dt);

        // Eqs. 21-23: tangential traction
        const scalar cf = sampler_.frictionCoefficient(i);
        const vector eTau = sampler_.tangentialDirection(i, U, h);
        const vector fTau = sampler_.tangentialForce(i, cf, eTau);

        // Per-node force on the body per unit density (area weight)
        vector nodeForce = (fN + fTau)*A[i];

        // Compressible: weight the node force by the local density at the
        // node's containing cell (mirrors the element multiplyForceRho)
        scalar rhoNode = 1.0;

        if (rhoPtr)
        {
            rhoNode = VGREAT;

            const label cellI = mesh_.findCell(X[i]);
            if (cellI >= 0)
            {
                rhoNode = (*rhoPtr)[cellI];
            }
            reduce(rhoNode, minOp<scalar>());

            if (!(rhoNode < VGREAT))
            {
                FatalErrorInFunction
                    << "Nacelle surface sample at " << X[i]
                    << " not found in mesh" << nl
                    << exit(FatalError);
            }

            nodeForce *= rhoNode;
        }

        // Per-unit-density distribution input (the kernel spread acts in the
        // global frame), and the SI on-body contract in the body frame
        sampler_.nodeForces_[i] = nodeForce;
        sampler_.forces_[i] =
            sampler_.rhoRef()*(sampler_.bodyToGlobal().T() & nodeForce);

        force_ += nodeForce;
        forceNormal_ += rhoNode*fN*A[i];
        forceTangential_ += rhoNode*fTau*A[i];
    }

    // Eq. 18: distribute the per-node forces onto the local cells. The
    // nacelle keeps the delivered exhaustive loop (nullptr candidate list)
    sampler_.distributeForce
    (
        forceField_,
        sampler_.nodeForces_,
        nodeH,
        nullptr
    );
}


void Foam::fv::nacelleSurfaceSource::addToEquation(fvMatrix<vector>& eqn)
{
    // Check dimensions of the force field and correct if necessary (mirrors
    // actuatorLineSource::addSup)
    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
    }

    eqn += forceField_;
}


void Foam::fv::nacelleSurfaceSource::writeOutput()
{
    if (writePerf_ and Pstream::master())
    {
        writePerf();

        if (writeNodePerf_)
        {
            writeNodePerf();
        }
    }
}


void Foam::fv::nacelleSurfaceSource::writePerf()
{
    const scalar rhoRef = sampler_.rhoRef();

    // Total force on the body in SI newtons
    const vector fSI = rhoRef*force_;
    const scalar fNMag = rhoRef*mag(forceNormal_);
    const scalar fTauMag = rhoRef*mag(forceTangential_);

    // cd = |F| / (0.5*rho*U^2*Aref)
    const scalar q = 0.5*rhoRef*sqr(sampler_.referenceVelocity())
                   * referenceArea_;
    const scalar cd = (q > VSMALL ? mag(fSI)/q : 0.0);

    *outputFile_ << mesh_.time().value()
        << "," << fSI.x() << "," << fSI.y() << "," << fSI.z()
        << "," << fNMag << "," << fTauMag << "," << cd << endl;
}


void Foam::fv::nacelleSurfaceSource::writeNodePerf()
{
    // Body-frame positions/normals and SI on-body forces (the contract)
    const List<point>& X = sampler_.positions();
    const List<vector>& n = sampler_.normals();
    const List<vector>& F = sampler_.forces();
    const List<scalar>& A = sampler_.areas();

    forAll(F, i)
    {
        *nodeFile_ << mesh_.time().value() << "," << i
            << "," << X[i].x() << "," << X[i].y() << "," << X[i].z()
            << "," << n[i].x() << "," << n[i].y() << "," << n[i].z()
            << "," << F[i].x() << "," << F[i].y() << "," << F[i].z()
            << "," << A[i] << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::nacelleSurfaceSource::nacelleSurfaceSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    sampler_(coeffs_, mesh),
    forceField_
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "force",
            dimForce/dimVolume/dimDensity,
            vector::zero
        )
    ),
    force_(vector::zero),
    forceNormal_(vector::zero),
    forceTangential_(vector::zero),
    referenceArea_(coeffs_.lookupOrDefault<scalar>("referenceArea", 1.0)),
    writePerf_(coeffs_.lookupOrDefault("writePerf", true)),
    writeNodePerf_(coeffs_.lookupOrDefault("writeNodePerf", false)),
    outputFile_(nullptr),
    nodeFile_(nullptr)
{
    read(dict_);

    if (writePerf_)
    {
        createOutputFile();
    }

    if (forceField_.writeOpt() == IOobject::AUTO_WRITE)
    {
        forceField_.write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::nacelleSurfaceSource::~nacelleSurfaceSource()
{
    delete outputFile_;
    delete nodeFile_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::nacelleSurfaceSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    calcForceField(eqn.psi());

    addToEquation(eqn);

    writeOutput();
}


void Foam::fv::nacelleSurfaceSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    calcForceField(eqn.psi(), &rho);

    addToEquation(eqn);

    writeOutput();
}


void Foam::fv::nacelleSurfaceSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // No-op: the nacelle actuator surface model injects no turbulence
    // source term (unlike the actuator line elements). Present for
    // interface completeness so axialFlowTurbineALSource::addSup compiles
    // and runs without effect.
}


bool Foam::fv::nacelleSurfaceSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        if (!coeffs_.readIfPresent("fieldNames", fieldNames_))
        {
            FatalErrorInFunction
                << "The nacelle surface source requires a 'fieldNames' entry"
                << nl << exit(FatalError);
        }
        applied_.setSize(fieldNames_.size(), false);

        referenceArea_ =
            coeffs_.lookupOrDefault<scalar>("referenceArea", 1.0);
        writePerf_ = coeffs_.lookupOrDefault("writePerf", true);
        writeNodePerf_ = coeffs_.lookupOrDefault("writeNodePerf", false);

        const bool writeForceField =
            coeffs_.lookupOrDefault("writeForceField", true);
        forceField_.writeOpt() =
            (writeForceField ? IOobject::AUTO_WRITE : IOobject::NO_WRITE);

        if (debug)
        {
            Info<< "Debugging for nacelleSurfaceSource on" << endl;
            printCoeffs();
        }

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::nacelleSurfaceSource::printCoeffs() const
{
    Info<< "Nacelle surface source properties:" << endl;
    Info<< "    nodes: " << sampler_.nNodes() << endl;
    Info<< "    reference velocity: " << sampler_.referenceVelocity() << endl;
    Info<< "    reference density: " << sampler_.rhoRef() << endl;
    Info<< "    reference area: " << referenceArea_ << endl;
}


// ************************************************************************* //
