/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 efirvida
-------------------------------------------------------------------------------
License
    This file is part of the solidBodyDisplacementLaplacianZone plugin for
    OpenFOAM.

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

#include "boundaryDecayDiffusivity.H"
#include "addToRunTimeSelectionTable.H"
#include "HashSet.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "patchWave.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundaryDecayDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        boundaryDecayDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::boundaryDecayDiffusivity::calcDecayFactor()
{
    const labelHashSet patchIDs
    (
        mesh().boundaryMesh().patchSet(patchNames_)
    );

    if (patchIDs.empty())
    {
        WarningInFunction
            << "No patches matched by " << patchNames_
            << ".  Decay factor set to 1 (no effect)." << endl;

        decayFactor_ = dimensionedScalar("one", dimless, 1.0);
        return;
    }

    // Compute cell-centre distance from selected patches via meshWave
    patchWave wave(mesh(), patchIDs, true);

    const scalarField& cellDist = wave.distance();

    // Build a cell-centred volScalarField for the decay factor
    volScalarField cellDecay
    (
        IOobject
        (
            "cellDecay",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        fvPatchFieldBase::zeroGradientType()
    );

    scalarField& cellDecayRef = cellDecay.primitiveFieldRef();

    forAll(cellDecayRef, celli)
    {
        const scalar xi = cellDist[celli] / decayDistance_;
        cellDecayRef[celli] = smoothStep(xi);
    }

    cellDecay.correctBoundaryConditions();

    // Interpolate to faces
    decayFactor_ = fvc::interpolate(cellDecay);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryDecayDiffusivity::boundaryDecayDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    motionDiffusivity(mesh),
    decayDistance_(readScalar(mdData)),
    patchNames_(mdData),
    basicDiffusivityPtr_(motionDiffusivity::New(mesh, mdData)),
    decayFactor_
    (
        IOobject
        (
            "decayFactor",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    )
{
    if (decayDistance_ <= SMALL)
    {
        FatalErrorInFunction
            << "decayDistance must be positive, got " << decayDistance_
            << exit(FatalError);
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryDecayDiffusivity::~boundaryDecayDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::boundaryDecayDiffusivity::operator()() const
{
    return decayFactor_ * basicDiffusivityPtr_->operator()();
}


void Foam::boundaryDecayDiffusivity::correct()
{
    basicDiffusivityPtr_->correct();
    calcDecayFactor();
}


// ************************************************************************* //
