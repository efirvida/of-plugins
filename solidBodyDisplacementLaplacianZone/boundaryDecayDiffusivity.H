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

Class
    Foam::boundaryDecayDiffusivity

Description
    Mesh motion diffusivity manipulator that multiplies a base diffusivity
    by a smooth decay factor near selected boundary patches.

    The decay factor transitions from 0 at the boundary patches to 1 at a
    user-specified distance, effectively preventing mesh deformation from
    reaching those patches.  This is particularly useful for AMI (Arbitrary
    Mesh Interface) boundaries where diffused deformation would break the
    interface connectivity.

    The decay factor is computed as a Hermite smooth-step:

    \f[
        f(\xi) = 3\xi^2 - 2\xi^3, \quad \xi = \min(d / d_{decay}, 1)
    \f]

    where \f$d\f$ is the cell distance to the nearest selected patch and
    \f$d_{decay}\f$ is the user-supplied decay distance.

Usage
    Example specification in \c dynamicMeshDict:

    \verbatim
        diffusivity  boundaryDecay 0.05 (AMI.*) quadratic inverseDistance (blade);
    \endverbatim

    where the entries are read sequentially from the Istream:

    \table
        Property     | Description                          | Required
        decayDistance | distance over which D ramps 0 -> 1  | yes
        patchNames   | patches to decay towards (wordRes)   | yes
        baseDiffusivity | any valid motionDiffusivity chain | yes
    \endtable

SourceFiles
    boundaryDecayDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryDecayDiffusivity_H
#define boundaryDecayDiffusivity_H

#include "motionDiffusivity.H"
#include "wordRes.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class boundaryDecayDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class boundaryDecayDiffusivity
:
    public motionDiffusivity
{
    // Private Data

        //- Distance over which the decay factor ramps from 0 to 1
        scalar decayDistance_;

        //- Patches to decay towards (supports regex)
        wordRes patchNames_;

        //- Wrapped base diffusivity
        autoPtr<motionDiffusivity> basicDiffusivityPtr_;

        //- Cached face-interpolated decay factor (0 at boundary, 1 far away)
        surfaceScalarField decayFactor_;


    // Private Member Functions

        //- Compute cell-centre distance from selected patches via meshWave
        //  and update the cached decayFactor_ field
        void calcDecayFactor();

        //- Hermite smooth-step: 3x^2 - 2x^3
        static inline scalar smoothStep(const scalar xi)
        {
            const scalar x = max(scalar(0), min(scalar(1), xi));
            return x * x * (scalar(3) - scalar(2) * x);
        }

        //- No copy construct
        boundaryDecayDiffusivity(const boundaryDecayDiffusivity&) = delete;

        //- No copy assignment
        void operator=(const boundaryDecayDiffusivity&) = delete;


public:

    //- Runtime type information
    TypeName("boundaryDecay");


    // Constructors

        //- Construct for the given fvMesh and data Istream
        boundaryDecayDiffusivity(const fvMesh& mesh, Istream& mdData);


    //- Destructor
    virtual ~boundaryDecayDiffusivity();


    // Member Functions

        //- Return diffusivity field (base * decayFactor)
        virtual tmp<surfaceScalarField> operator()() const;

        //- Correct the motion diffusivity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
