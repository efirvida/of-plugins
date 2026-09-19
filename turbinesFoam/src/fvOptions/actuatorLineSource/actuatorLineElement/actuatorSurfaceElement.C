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

#include "actuatorSurfaceElement.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"
#include "pointField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorSurfaceElement, 0);
    addToRunTimeSelectionTable
    (
        actuatorLineElement,
        actuatorSurfaceElement,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::vector Foam::fv::actuatorSurfaceElement::chordPoint(label k) const
{
    // Midpoint of chord strip k: f_k = (k + 0.5)/nChordwise_
    scalar fk = (k + 0.5)/nChordwise_;
    vector chordUnit = chordDirection_/mag(chordDirection_);

    return position_ + (chordMount_ - fk)*chordLength_*chordUnit;
}


Foam::scalar Foam::fv::actuatorSurfaceElement::calcProjectionEpsilon()
{
    // Lookup Gaussian coeffs from profileData dict if present
    dictionary GaussianCoeffs = profileData_.dict().subOrEmptyDict
    (
        "GaussianCoeffs"
    );
    scalar meshFactor = GaussianCoeffs.lookupOrDefault("meshFactor", 2.0);

    // Projection width based on local cell size only (from Troldborg (2008)).
    // The actuator surface model couples to fine meshes, so no chord-length
    // term enters the epsilon calculation.
    scalar epsilon = VGREAT;
    const scalarField& V = mesh_.V();
    label posCellI = findCell(position_);
    if (posCellI >= 0)
    {
        epsilon = 2.0*Foam::cbrt(V[posCellI]);
        epsilon *= meshFactor; // Cell could have non-unity aspect ratio
    }

    // Reduce epsilon over all processors
    reduce(epsilon, minOp<scalar>());

    // If epsilon is not reduced, position is not in the mesh
    if (not (epsilon < VGREAT))
    {
        // Raise fatal error since mesh size cannot be detected
        FatalErrorIn("void actuatorSurfaceElement::calcProjectionEpsilon()")
            << "Position of " << name_ << " not found in mesh"
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "    epsilon (mesh-based): " << epsilon << endl;
    }

    return epsilon;
}


void Foam::fv::actuatorSurfaceElement::applyForceField
(
    volVectorField& forceField
)
{
    // Calculate projection width
    scalar epsilon = calcProjectionEpsilon();
    scalar projectionRadius = (epsilon*Foam::sqrt(Foam::log(1.0/0.001)));

    // Each chord strip carries an equal share of the element's total force
    vector stripForce = forceVector_/nChordwise_;

    // Bounding box around the chord line, inflated by the projection radius,
    // used as a prefilter for the cell loop
    vector chordUnit = chordDirection_/mag(chordDirection_);
    pointField chordEndpoints(2);
    chordEndpoints[0] = position_ + chordMount_*chordLength_*chordUnit;
    chordEndpoints[1] = position_ - (1.0 - chordMount_)*chordLength_*chordUnit;
    boundBox chordBox(chordEndpoints, false);
    chordBox.inflate(projectionRadius);

    if (debug)
    {
        Info<< "    nChordwise: " << nChordwise_ << endl;
        Info<< "    stripForce: " << stripForce << endl;
        Info<< "    sphereRadius: " << chordLength_ + projectionRadius << endl;
    }

    // Apply force to the cells within the chord line's sphere of influence
    forAll(mesh_.cells(), cellI)
    {
        const point& C = mesh_.C()[cellI];
        if (chordBox.contains(C))
        {
            for (label k = 0; k < nChordwise_; k++)
            {
                scalar dis = mag(C - chordPoint(k));
                if (dis <= chordLength_ + projectionRadius)
                {
                    scalar factor = Foam::exp(-Foam::sqr(dis/epsilon))
                                  / (Foam::pow(epsilon, 3)
                                  * Foam::pow(Foam::constant::mathematical::pi, 1.5));
                    // forceField is opposite forceVector
                    forceField[cellI] += -stripForce*factor;
                }
            }
        }
    }
}


void Foam::fv::actuatorSurfaceElement::calculateInflowVelocity
(
    const volVectorField& Uin
)
{
    // Find local flow velocity by interpolating to the midpoint of each
    // chord strip, then average over all strips
    inflowVelocity_ = vector(VGREAT, VGREAT, VGREAT);
    vector velocitySum = vector(0.0, 0.0, 0.0);
    interpolationCellPoint<vector> UInterp(Uin);

    for (label k = 0; k < nChordwise_; k++)
    {
        vector sampleVelocity = vector(VGREAT, VGREAT, VGREAT);
        vector samplePoint = chordPoint(k);

        // Sample the velocity
        label sampleCellI = findCell(samplePoint);
        if (sampleCellI >= 0)
        {
            sampleVelocity = UInterp.interpolate
            (
                samplePoint,
                sampleCellI
            );
        }

        // Reduce inflow velocity over all processors
        reduce(sampleVelocity, minOp<vector>());

        // If inflow velocity is not detected, position is not in the mesh
        if (not (sampleVelocity[0] < VGREAT))
        {
            // Raise fatal error since inflow velocity cannot be detected
            FatalErrorIn("void actuatorSurfaceElement::calculateInflowVelocity()")
                << "Inflow velocity point for " << name_
                << " not found in mesh"
                << abort(FatalError);
        }

        velocitySum = velocitySum + sampleVelocity;
    }

    // Set inflow velocity as the mean over the chord strips
    inflowVelocity_ = 1.0/nChordwise_*velocitySum;

    if (debug)
    {
        Info<< "    chord-averaged inflow velocity: " << inflowVelocity_
            << " (" << nChordwise_ << " chordwise strips)" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorSurfaceElement::actuatorSurfaceElement
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    actuatorLineElement(name, dict, mesh),
    nChordwise_(dict.lookupOrDefault<label>("nChordwise", 5))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorSurfaceElement::~actuatorSurfaceElement()
{}


// ************************************************************************* //