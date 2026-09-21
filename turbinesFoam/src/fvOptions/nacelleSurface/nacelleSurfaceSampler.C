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
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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
    surfaceSamplerBase(dict, mesh),
    forces_(),
    nodeForces_(),
    referenceVelocity_(dict.get<scalar>("referenceVelocity")),
    nu_(-1.0),
    cfOverride_(dict.lookupOrDefault<scalar>("cf", -1.0)),
    streamwiseDirection_(vector(1, 0, 0)),
    noseStreamwise_(0.0),
    streamwiseCoords_()
{
    forces_.setSize(nNodes(), vector::zero);
    nodeForces_.setSize(nNodes(), vector::zero);

    // The nacelle axis is the streamwise direction (the nacelle is aligned
    // with the incoming flow by construction). createBodyFrame() has already
    // validated the entry and completed the orthonormal basis; repeat the
    // exact normalization to keep the streamwise coordinates identical
    streamwiseDirection_ = dict.lookupOrDefault("bodyAxis", vector(1, 0, 0));
    streamwiseDirection_ /= mag(streamwiseDirection_);

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
    streamwiseCoords_.setSize(nNodes());
    noseStreamwise_ = VGREAT;

    forAll(positions_, i)
    {
        streamwiseCoords_[i] = positions_[i] & streamwiseDirection_;
        noseStreamwise_ = min(noseStreamwise_, streamwiseCoords_[i]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::nacelleSurfaceSampler::~nacelleSurfaceSampler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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


// ************************************************************************* //
