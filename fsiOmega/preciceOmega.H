/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2026 OpenFOAM Foundation
    Copyright (C) 2024-2026 Contributors to the fsiOmega library
-------------------------------------------------------------------------------
License
    This file is part of the fsiOmega library for OpenFOAM.

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
    Foam::Function1Types::preciceOmega

Group
    grpFsiFunction1

Description
    Angular velocity function that reads omega from a uniformDimensionedScalarField
    updated by the preCICE adapter. This enables FSI coupling where the angular
    velocity is computed by an external structural solver and communicated via preCICE.

    \heading How It Works

    The preCICE OpenFOAM adapter receives angular velocity data from a coupled
    participant (e.g., a structural solver) and stores it in a 
    uniformDimensionedScalarField. This Function1 reads that field value
    and provides it to the solidBodyMotion solver.

    \heading Data Flow

    \verbatim
    External Solver (e.g., MBDyn, structural FEM)
        |
        v
    preCICE (coupling library)
        |
        v
    OpenFOAM Adapter (AngularVelocity reader)
        |
        v
    uniformDimensionedScalarField "omega"
        |
        v
    preciceOmega Function1
        |
        v
    solidBodyMotion (rotatingMotion)
    \endverbatim

    \heading Usage

    In constant/dynamicMeshDict:
    \verbatim
    dynamicFvMesh       dynamicMotionSolverFvMesh;
    motionSolverLibs    ("libsolidBodyDisplacementLaplacianZoneFvMotionSolver.so" "libfvMotionSolvers.so" "libfsiOmega.so");
    motionSolver        solidBodyDisplacementLaplacianZone;

    solidBodyDisplacementLaplacianZoneCoeffs
    {
        cellZone    rotorZone;

        solidBodyMotionFunction rotatingMotion;

        rotatingMotionCoeffs
        {
            origin      (0 0 0);
            axis        (0 0 1);

            omega       preciceOmega;
            preciceOmegaCoeffs
            {
                // Name of the field updated by preCICE adapter
                // Must match nameOmegaField in preciceDict FSI section
                fieldName   omega;
            }
        }

        diffusivity quadratic inverseDistance (rotorTip);
    }
    \endverbatim

    \heading Configuration in preciceDict

    The preCICE adapter must be configured to read AngularVelocity:
    \verbatim
    FSI
    {
        nameOmegaField  omega;      // Must match fieldName above
        rotationAxis    (0 0 1);
    }

    interfaces
    {
        Interface1
        {
            mesh        Fluid-Mesh;
            patches     (rotor);
            locations   faceCenters;
            
            readData    (AngularVelocity);
            writeData   (Force);
        }
    };
    \endverbatim

See also
    Foam::Function1
    preciceAdapter::FSI::AngularVelocity

SourceFiles
    preciceOmega.C

\*---------------------------------------------------------------------------*/

#ifndef Function1Types_preciceOmega_H
#define Function1Types_preciceOmega_H

#include "Function1.H"
#include "objectRegistry.H"
#include "uniformDimensionedFields.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                         Class preciceOmega Declaration
\*---------------------------------------------------------------------------*/

class preciceOmega
:
    public Function1<scalar>
{
    // Private Data

        //- Reference to the object registry (mesh)
        const objectRegistry* obrPtr_;

        //- Name of the uniformDimensionedScalarField to read from
        word fieldName_;

        //- Angular velocity field (non-owning pointer)
        mutable uniformDimensionedScalarField* omegaField_;

        //- Angular velocity field (owning pointer, we bind to omegaField_)
        mutable std::unique_ptr<uniformDimensionedScalarField> omegaFieldOwning_;


    // Private Member Functions

        //- Read parameters from dictionary
        void read(const dictionary& dict);

        //- Initialize the omega field (find existing or create new)
        void initField() const;

        //- No copy assignment
        void operator=(const preciceOmega&) = delete;


public:

    //- Runtime type information
    TypeName("preciceOmega");


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        preciceOmega
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy constructor
        explicit preciceOmega(const preciceOmega& rhs);

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const
        {
            return tmp<Function1<scalar>>(new preciceOmega(*this));
        }


    //- Destructor
    virtual ~preciceOmega() = default;


    // Member Functions

        //- Return the angular velocity value from the preCICE-updated field
        virtual scalar value(const scalar t) const;

        //- Integrate between two time values
        //  Note: Since omega can change at each coupling step,
        //  integration is approximate using the current value
        virtual scalar integrate(const scalar t1, const scalar t2) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
