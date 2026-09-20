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

#include "nacelleSurfaceSampler.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "dimensionedScalar.H"
#include "mathematicalConstants.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::nacelleSurfaceSampler::createBodyFrame
(
    const dictionary& dict
)
{
    const vector axis = dict.lookupOrDefault("bodyAxis", vector(1, 0, 0));

    if (mag(axis) < SMALL)
    {
        FatalErrorInFunction
            << "The nacelle 'bodyAxis' entry must be non-zero" << nl
            << exit(FatalError);
    }

    // Body axis e1 is the nacelle axis expressed in the global frame; complete
    // the orthonormal basis deterministically
    const vector e1 = axis/mag(axis);

    vector e2 = e1 ^ vector(0, 0, 1);
    if (mag(e2) < SMALL)
    {
        e2 = e1 ^ vector(0, 1, 0);
    }
    e2 /= mag(e2);

    const vector e3 = e1 ^ e2;

    // Columns of bodyToGlobal_ are the body axes expressed in the global frame
    bodyToGlobal_ = tensor
    (
        e1.x(), e2.x(), e3.x(),
        e1.y(), e2.y(), e3.y(),
        e1.z(), e2.z(), e3.z()
    );

    // The nacelle axis is the streamwise direction (the nacelle is aligned
    // with the incoming flow by construction)
    streamwiseDirection_ = e1;

    identityBodyFrame_ = (bodyOrigin_ == vector::zero)
                      && (e1 == vector(1, 0, 0));

    if (identityBodyFrame_)
    {
        positionsBody_ = positions_;
        normalsBody_ = normals_;
    }
    else
    {
        const tensor globalToBody = bodyToGlobal_.T();

        forAll(positions_, i)
        {
            positionsBody_[i] = globalToBody & (positions_[i] - bodyOrigin_);
            normalsBody_[i] = globalToBody & normals_[i];
        }
    }
}


void Foam::fv::nacelleSurfaceSampler::readViscosity(const dictionary& dict)
{
    nu_ = dict.lookupOrDefault<scalar>("nu", -1.0);

    if (nu_ >= 0.0)
    {
        return;
    }

    // Fall back to constant/transportProperties (the openfoam convention)
    const fileName tpPath = mesh_.time().constant()/"transportProperties";

    if (isFile(tpPath))
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        );

        if (transportProperties.found("nu"))
        {
            nu_ = dimensionedScalar
            (
                "nu",
                dimViscosity,
                transportProperties
            ).value();
        }
    }

    if (nu_ <= 0.0)
    {
        FatalErrorInFunction
            << "The nacelle surface source requires a positive kinematic "
            << "viscosity: set 'nu' in the nacelleSurfaceSourceCoeffs or "
            << "provide a 'nu' entry in constant/transportProperties" << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::nacelleSurfaceSampler::nacelleSurfaceSampler
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    surface_(),
    positions_(),
    normals_(),
    areas_(),
    positionsBody_(),
    normalsBody_(),
    forces_(),
    nodeForces_(),
    bodyOrigin_(dict.lookupOrDefault("bodyOrigin", vector::zero)),
    bodyToGlobal_(tensor::I),
    identityBodyFrame_(true),
    referenceVelocity_(dict.get<scalar>("referenceVelocity")),
    nu_(-1.0),
    cfOverride_(dict.lookupOrDefault<scalar>("cf", -1.0)),
    rhoRef_(dict.lookupOrDefault<scalar>("rho", 1.0)),
    streamwiseDirection_(vector(1, 0, 0)),
    noseStreamwise_(0.0),
    streamwiseCoords_()
{
    if (!dict.found("geometry"))
    {
        FatalErrorInFunction
            << "The nacelle surface source requires a 'geometry' entry "
            << "giving the surface triangulation file" << nl
            << exit(FatalError);
    }

    fileName geometryPath = dict.get<fileName>("geometry");

    if (!isFile(geometryPath))
    {
        // Also accept a path relative to the case directory
        const fileName casePath = mesh_.time().path()/geometryPath;

        if (isFile(casePath))
        {
            geometryPath = casePath;
        }
    }

    if (!isFile(geometryPath))
    {
        FatalErrorInFunction
            << "Nacelle surface file " << geometryPath << " not found" << nl
            << exit(FatalError);
    }

    // triSurface::New auto-detects the ASCII/binary STL format
    surface_.reset(triSurface::New(geometryPath));

    if (surface_->size() == 0)
    {
        FatalErrorInFunction
            << "Nacelle surface file " << geometryPath
            << " contains no triangles" << nl
            << exit(FatalError);
    }

    // One node per triangle: centroid position, outward unit normal, area
    positions_ = surface_->faceCentres();
    normals_ = surface_->faceNormals();
    areas_ = surface_->magFaceAreas();

    positionsBody_.setSize(positions_.size());
    normalsBody_.setSize(normals_.size());
    forces_.setSize(positions_.size(), vector::zero);
    nodeForces_.setSize(positions_.size(), vector::zero);

    createBodyFrame(dict);
    readViscosity(dict);

    // Validate the friction model selection: 'constant' requires a cf value
    const word cfModel = dict.lookupOrDefault<word>("cfModel", "schultzGrunow");

    if (cfModel == "constant")
    {
        if (cfOverride_ < 0.0)
        {
            FatalErrorInFunction
                << "cfModel 'constant' requires a non-negative 'cf' entry"
                << nl << exit(FatalError);
        }
    }
    else if (cfModel != "schultzGrunow")
    {
        FatalErrorInFunction
            << "Unknown cfModel '" << cfModel << "'" << nl
            << "Valid models: schultzGrunow, constant" << nl
            << exit(FatalError);
    }

    // Streamwise coordinate of every node and of the most-upstream node (the
    // nose), used by the Schultz-Grunow friction relation
    streamwiseCoords_.setSize(positions_.size());
    noseStreamwise_ = VGREAT;

    forAll(positions_, i)
    {
        streamwiseCoords_[i] = positions_[i] & streamwiseDirection_;
        noseStreamwise_ = min(noseStreamwise_, streamwiseCoords_[i]);
    }

    Info<< "Nacelle surface sampler: read " << positions_.size()
        << " triangles from " << geometryPath << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::nacelleSurfaceSampler::~nacelleSurfaceSampler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::nacelleSurfaceSampler::cellSize(const point& X) const
{
    scalar h = VGREAT;

    const label cellI = mesh_.findCell(X);

    if (cellI >= 0)
    {
        h = Foam::cbrt(mesh_.V()[cellI]);
    }

    // Reduce the sentinel over all processors; the rank owning the containing
    // cell contributes its local cell size
    reduce(h, minOp<scalar>());

    if (!(h < VGREAT))
    {
        FatalErrorInFunction
            << "Nacelle surface sample at " << X << " not found in mesh" << nl
            << exit(FatalError);
    }

    return h;
}


Foam::vector Foam::fv::nacelleSurfaceSampler::normalForce
(
    const label i,
    const vector& uTilde,
    const scalar h,
    const scalar dt
) const
{
    const vector& eN = normals_[i];

    // Eq. 19 with the desired surface velocity u^d = 0 (stationary nacelle).
    // With the paper's signs as written this is the force of the flow on the
    // surface (force on the body); the source applies -F_i to the flow.
    return (h*(uTilde & eN)/dt)*eN;
}


Foam::vector Foam::fv::nacelleSurfaceSampler::tangentialForce
(
    const label,
    const scalar cf,
    const vector& eTau
) const
{
    // Eq. 21; the node index is part of the per-node interface, the reference
    // velocity is a global value
    return 0.5*cf*sqr(referenceVelocity_)*eTau;
}


Foam::scalar Foam::fv::nacelleSurfaceSampler::frictionCoefficient
(
    const label i
) const
{
    if (cfOverride_ >= 0.0)
    {
        return cfOverride_;
    }

    // Eq. 22: Schultz-Grunow, keyed by the streamwise distance behind the nose
    const scalar dist = max(streamwiseCoords_[i] - noseStreamwise_, 0.0);
    const scalar Rex = referenceVelocity_*dist/nu_;

    // The relation is only defined for log(Rex) > 0. In the nose region
    // Rex -> 0 and the zero-pressure-gradient assumption is invalid
    // (documented); return zero friction there instead of a complex value
    if (Rex <= 1.0)
    {
        return 0.0;
    }

    return 0.37*Foam::pow(Foam::log(Rex), -2.584);
}


Foam::vector Foam::fv::nacelleSurfaceSampler::tangentialDirection
(
    const label i,
    const volVectorField& U,
    const scalar h
) const
{
    // Eq. 23: normalized velocity at the off-wall probe point X + h*e_n
    const point probe = positions_[i] + h*normals_[i];

    const vector uProbe = interpolateVelocity(probe, U, h);
    const scalar magU = mag(uProbe);

    if (magU < SMALL)
    {
        // Stagnation-point guard: e_tau is undefined for zero probe velocity
        return vector::zero;
    }

    return uProbe/magU;
}


Foam::scalar Foam::fv::nacelleSurfaceSampler::kernel(const scalar r)
{
    // Eq. 8: smoothed four-point cosine kernel, support |r| <= 2.5
    const scalar ar = mag(r);

    if (ar <= 1.5)
    {
        return
            0.25
          + Foam::sin
            (
                constant::mathematical::pi*(2.0*ar + 1.0)/4.0
            )/(2.0*constant::mathematical::pi)
          - Foam::sin
            (
                constant::mathematical::pi*(2.0*ar - 1.0)/4.0
            )/(2.0*constant::mathematical::pi);
    }
    else if (ar <= 2.5)
    {
        return
            0.625
          - 0.25*ar
          - Foam::sin
            (
                constant::mathematical::pi*(2.0*ar - 1.0)/4.0
            )/(2.0*constant::mathematical::pi);
    }

    return 0.0;
}


Foam::vector Foam::fv::nacelleSurfaceSampler::interpolateVelocity
(
    const point& X,
    const volVectorField& U,
    const scalar h
) const
{
    // Eq. 7: kernel sum over the local cells within the kernel support,
    // delta_h*V = phi_x*phi_y*phi_z for a uniform local cell size h
    vector sumU = vector::zero;
    scalar sumW = 0.0;
    const scalar radius = 2.5*h;

    forAll(mesh_.cells(), cellI)
    {
        const vector d = mesh_.C()[cellI] - X;

        // Bounding-box prefilter on the kernel support
        if
        (
            mag(d.x()) > radius
         || mag(d.y()) > radius
         || mag(d.z()) > radius
        )
        {
            continue;
        }

        const scalar w =
            kernel(d.x()/h)*kernel(d.y()/h)*kernel(d.z()/h);

        if (w > 0.0)
        {
            sumU += w*U[cellI];
            sumW += w;
        }
    }

    // Every rank holds the full node list but only its local cells, so the
    // kernel sums are reduced globally
    returnReduce(sumW, sumOp<scalar>());
    returnReduce(sumU, sumOp<vector>());

    if (sumW < SMALL)
    {
        // No cells within the kernel support (e.g. an off-wall probe outside
        // the mesh); treated as a zero velocity
        return vector::zero;
    }

    return sumU/sumW;
}


// ************************************************************************* //
