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
    Foam::fv::bladeSurfaceSource

Description
    Per-blade distributor of the BEM element forces over an imported blade
    triangulation (blade actuator surface model, S2; design D1/D8/D9).

    A plain helper owned by the blade \c actuatorLineSource when the blade
    subdictionary carries \c surfaceGeometry; it is not an \c fv::option and
    is not added to the option run-time selection table (D1). It owns a
    \c bladeSurfaceSampler and consumes the owning blade's
    \c PtrList<actuatorLineElement> after the element loop:

    - per node: \c F_node = patchAreaShare_i*elements_[patch_i].force()
      (per unit density, force on the blade)
    - per candidate cell: \c forceField_[c] -= F_node*(w/V_c), the paper's
      Eq. 18 with the S1 sign, over the node's bounded candidate list
    - the returned total is the distributed force (global, \c returnReduce);
      partition of unity makes it equal the summed element forces in the
      incompressible reference, and it replaces the element-loop total so the
      reported blade force matches the applied load (D8-consistent)

    Unit convention: the element force is per unit density. The SI contract
    (sampler \c forces(), per-station CSV) uses \c rho (default 1) as the
    reference density, the S1 convention. \c moment() returns the moment of
    the distributed load about a point, per unit density, global frame.

    Compressible overload: the element's public \c force() already carries the
    element-local density (\c multiplyForceRho). The distributor recovers the
    per-unit-density share and re-applies the density at the node, following
    the S1 node-local density convention; the incompressible reference
    (Phase VI) has a unit ratio.

    Output (master rank only, \c postProcessing/bladeSurface/):

    \verbatim
    <owner>.surface.csv         per-station (writePerf, default true)
    <owner>.surface_nodes.csv   per-node (writeNodePerf, default false)
    \endverbatim

    The per-station schema is
    \c time,station,root_dist,area,force_x,force_y,force_z,c_ref_n,c_ref_t,f_ref_n,f_ref_t
    and the per-node schema is
    \c time,node,x,y,z,nx,ny,nz,fx,fy,fz,area,station,chord_fraction
    (design section 4.5); the per-node quantities are the body-frame SI
    contract.

SourceFiles
    bladeSurfaceSource.C

\*---------------------------------------------------------------------------*/

#ifndef bladeSurfaceSource_H
#define bladeSurfaceSource_H

#include "bladeSurfaceSampler.H"
#include "volFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

class actuatorLineElement;

/*---------------------------------------------------------------------------*\
                  Class bladeSurfaceSource Declaration
\*---------------------------------------------------------------------------*/

class bladeSurfaceSource
{
    // Private data

        //- Surface sampler (geometry, partition, candidates, node forces)
        bladeSurfaceSampler sampler_;

        //- Mesh the surface is sampled on
        const fvMesh& mesh_;

        //- Owning blade elements (the BEM force source)
        PtrList<actuatorLineElement>& elements_;

        //- Number of element patches
        label nElements_;

        //- Output file base name (<ownerName>.surface)
        word name_;

        //- Switch for the per-station CSV
        bool writePerf_;

        //- Switch for the opt-in per-node CSV
        bool writeNodePerf_;

        //- Per-station output file stream
        OFstream* stationFile_;

        //- Per-node output file stream
        OFstream* nodeFile_;


    // Private Member Functions

        //- No copy construct
        bladeSurfaceSource(const bladeSurfaceSource&) = delete;

        //- No copy assignment
        void operator=(const bladeSurfaceSource&) = delete;

        //- Create the output file(s) under postProcessing/bladeSurface
        void createOutputFiles();

        //- Write the CSV(s) on the master rank
        void writeOutput();

        //- Write the per-station CSV
        void writeStationCsv();

        //- Write the per-node CSV
        void writeNodeCsv();


public:

    // Constructors

        //- Construct from the owning blade name, the blade subdictionary,
        //  the mesh and the owning blade elements
        bladeSurfaceSource
        (
            const word& ownerName,
            const dictionary& dict,
            const fvMesh& mesh,
            PtrList<actuatorLineElement>& elements
        );


    //- Destructor
    ~bladeSurfaceSource();


    // Member Functions

        // Distribution

            //- Distribute the element forces over the surface nodes (Eq. 18)
            //  into the owning blade's force field. \c rhoPtr selects the
            //  compressible overload (local-density scaling); \c bladeForce
            //  receives the distributed total on return.
            void distribute
            (
                volVectorField& bladeForceField,
                vector& bladeForce,
                const volScalarField* rhoPtr = nullptr
            );


        // Evaluation

            //- Moment of the distributed load about point, per unit density,
            //  global frame (D8; replaces the element moment when active)
            vector moment(const point& p) const;


        // Access

            //- Return the surface sampler
            inline const bladeSurfaceSampler& sampler() const;

            //- A constructed distributor is active: the owning blade source
            //  controls activation by construction (D1)
            inline bool active() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bladeSurfaceSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
