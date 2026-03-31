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
    Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Single motion solver that combines rigid solid-body rotation of a selected
    mesh zone (cellZone or cellSet) with FSI deformation propagation via
    Laplacian diffusion.

    This solver implements the FSI motion model:

    \f[
        U_{total}(x) = R(t) \cdot (x + u_{fsi}(x)) - x
    \f]

    where:
    - \f$x\f$ are the original point coordinates
    - \f$u_{fsi}\f$ is the structural deformation (applied on the patch specified
      in the diffusivity definition, propagated diffusively within the rotation zone)
    - \f$R(t)\f$ is the rigid rotation transformation from solidBodyMotionFunction

    The deformation and rotation are applied in the correct physical order:
    1. Apply structural deformation to the diffusivity reference patch
    2. Propagate deformation diffusively within the rotation zone
    3. Apply rigid rotation to the entire deformed zone

    For overset meshes, the deformation stays confined to the rotation zone
    and does NOT propagate to acceptor cells.

    \heading Usage in constant/dynamicMeshDict

    \verbatim
    dynamicFvMesh   dynamicOversetZoneDisplacementFvMesh;  // for overset
    // dynamicFvMesh   dynamicMotionSolverFvMesh;           // without overset

    motionSolverLibs
    (
        "libsolidBodyDisplacementLaplacianZoneFvMotionSolver.so"
        "libfvMotionSolvers.so"
    );

    motionSolver    solidBodyDisplacementLaplacianZone;

    solidBodyDisplacementLaplacianZoneCoeffs
    {
        cellZone        rotorZone;      // Zone that rotates as a rigid body

        solidBodyMotionFunction  rotatingMotion;
        rotatingMotionCoeffs
        {
            origin      (0 0 0);
            axis        (0 1 0);
            omega       158;            // rad/s, or Function1<scalar>
        }

        // Diffusivity: the patch name in parentheses is the reference patch
        // where FSI deformation (u_fsi) is applied as fixedValue BC
        diffusivity quadratic inverseDistance (bladeSurface);

        // Optional: fix specific points
        // frozenPointsZone frozenZone;
    }
    \endverbatim

    \heading FSI / preCICE usage

    The structural deformation \f$u_{fsi}\f$ is set on the cellDisplacement field
    as a fixedValue boundary condition on the diffusivity reference patch.
    The preCICE adapter (or any other mechanism) updates this field before
    each mesh update.

    Angular velocity coupling via preCICE:

    \verbatim
        omega       preciceOmega;
        preciceOmegaCoeffs
        {
            fieldName   omega;   // must match nameOmegaField in preciceDict
        }
    \endverbatim

See also
    Foam::solidBodyDisplacementLaplacianFvMotionSolver
    Foam::displacementLaplacianFvMotionSolver
    Foam::dynamicOversetZoneDisplacementFvMesh

SourceFiles
    solidBodyDisplacementLaplacianZoneFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyDisplacementLaplacianZoneFvMotionSolver_H
#define solidBodyDisplacementLaplacianZoneFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class motionInterpolation;
class motionDiffusivity;
class solidBodyMotionFunction;

/*---------------------------------------------------------------------------*\
     Class solidBodyDisplacementLaplacianZoneFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class solidBodyDisplacementLaplacianZoneFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Motion function for rigid transformation
        autoPtr<solidBodyMotionFunction> solidBodyMotionPtr_;

        //- Cell-centre motion field (receives u_fsi from FSI solver via BC)
        mutable volVectorField cellDisplacement_;

        //- Optional point-position field (position boundary conditions)
        mutable autoPtr<pointVectorField> pointLocation_;

        //- Interpolation from cell displacement to points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Diffusivity controlling motion propagation
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Index of the frozen-points zone (-1 = none)
        label frozenPointsZone_;

        //- Indices of cells belonging to the selected zone
        labelList cellIDs_;

        //- Reference cell centres for absolute rigid displacement evaluation
        vectorField cellCentres0_;

        //- Indices of points belonging to the selected zone
        labelList pointIDs_;

        //- True when rigid transform is applied to the whole mesh
        bool moveAllCells_;

        //- Threshold for identifying overset hole cells (cellMask < this)
        static constexpr scalar holeCellThreshold_ = 0.5;


    // Private Member Functions

        //- No copy construct
        solidBodyDisplacementLaplacianZoneFvMotionSolver
        (
            const solidBodyDisplacementLaplacianZoneFvMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const solidBodyDisplacementLaplacianZoneFvMotionSolver&
        ) = delete;

        //- Common initialization: resolve zone/set, collect point IDs
        void initZoneAndPointIDs();


public:

    //- Runtime type information
    TypeName("solidBodyDisplacementLaplacianZone");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        solidBodyDisplacementLaplacianZoneFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        solidBodyDisplacementLaplacianZoneFvMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    ~solidBodyDisplacementLaplacianZoneFvMotionSolver();


    // Member Functions

        //- Return reference to the diffusivity field
        motionDiffusivity& diffusivity();

        //- Return point location from current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
