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

Class
    Foam::fv::nacelleSurfaceSource

Description
    Cell-set momentum source implementing the nacelle/hub actuator surface
    model of Yang & Sotiropoulos (arXiv:1702.02108v4, Sec. 2.2).

    The source owns a \c nacelleSurfaceSampler, which reads the triangulated
    nacelle surface from STL and computes per-node normal (Eq. 19) and
    tangential (Eqs. 21-23) tractions per unit density. Each node force is
    area-weighted and distributed onto the background cells with the paper's
    smoothed four-point cosine kernel (Eq. 8, Eq. 18), scaled by \c 1/V_cell,
    which preserves the total force by the kernel's partition of unity.

    Usage (standalone, rotor-less validation case):

    \verbatim
    nacelle
    {
        type            nacelleSurfaceSource;
        active          true;

        nacelleSurfaceSourceCoeffs
        {
            fieldNames          (U);
            selectionMode       cellSet;
            cellSet             nacelleCells;

            geometry            "geometry/nacelle.stl";
            referenceVelocity   1.0;
            rho                 1.0;
            nu                  1e-3;
            cfModel             schultzGrunow;  // or constant
            cf                  -1.0;
            bodyOrigin          (0 0 0);
            bodyAxis            (1 0 0);
            referenceArea       3.14159265;
            writeForceField     true;
            writePerf           true;
            writeNodePerf       false;
        }
    }
    \endverbatim

    The source is also composed by
    \c axialFlowTurbineALSource::createNacelle() as an owned child
    (\c nacelle_ ), constructed from the case's \c nacelle {} subdictionary.

    Unit convention (paper and \c cellSetOption ): the incompressible source
    is force per unit volume per unit density (m/s²). The compressible
    overload multiplies the per-node force and the field by the local
    density. \c forces() of the sampler reports the SI force on the body.

SourceFiles
    nacelleSurfaceSource.C

\*---------------------------------------------------------------------------*/

#ifndef nacelleSurfaceSource_H
#define nacelleSurfaceSource_H

#include "cellSetOption.H"
#include "volFields.H"
#include "OFstream.H"
#include "nacelleSurfaceSampler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class nacelleSurfaceSource Declaration
\*---------------------------------------------------------------------------*/

class nacelleSurfaceSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Surface sampler (geometry, force model, kernel)
        nacelleSurfaceSampler sampler_;

        //- Force field from all nodes (per unit density, "force." + name_)
        volVectorField forceField_;

        //- Total force on the body per unit density
        vector force_;

        //- Normal contribution to the total force per unit density
        vector forceNormal_;

        //- Tangential contribution to the total force per unit density
        vector forceTangential_;

        //- Reference area for the drag coefficient output
        scalar referenceArea_;

        //- Switch for writing the total force CSV
        bool writePerf_;

        //- Switch for writing the per-node CSV
        bool writeNodePerf_;

        //- Total force output file stream
        OFstream* outputFile_;

        //- Per-node output file stream
        OFstream* nodeFile_;


    // Protected Member Functions

        //- Create the performance output file(s)
        virtual void createOutputFile();

        //- Compute the per-node forces and fill forceField_ (Eqs. 18, 19,
        //  21-23); rho is used for the compressible overload
        void calcForceField
        (
            const volVectorField& U,
            const volScalarField* rhoPtr = nullptr
        );

        //- Distribute the per-node forces onto the local cells (Eq. 18)
        void distributeForce
        (
            volVectorField& ff,
            const List<scalar>& nodeH
        );

        //- Add the force field to the momentum equation
        void addToEquation(fvMatrix<vector>& eqn);

        //- Write the output CSV(s) on the master rank
        void writeOutput();

        //- Write the total-force CSV
        void writePerf();

        //- Write the per-node CSV
        void writeNodePerf();


public:

    //- Runtime type information
    TypeName("nacelleSurfaceSource");


    // Constructors

        //- Construct from components
        nacelleSurfaceSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~nacelleSurfaceSource();


    // Member Functions

        // Access

            //- Return the surface sampler
            inline const nacelleSurfaceSampler& sampler() const;

            //- Return the total force on the body per unit density
            inline const vector& force() const;

            //- Return the force field
            inline const volVectorField& forceField() const;


        // Source term addition

            //- Add source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add source term to turbulence scalars (no-op: the nacelle
            //  model injects no turbulence, unlike the actuator line elements)
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);

            //- Print dictionary values
            virtual void printCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nacelleSurfaceSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
