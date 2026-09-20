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
    Foam::fv::nacelleSurfaceSampler

Description
    Reusable surface-sampling object for the nacelle/hub actuator surface
    model (Yang & Sotiropoulos, arXiv:1702.02108v4, Sec. 2.2).

    The nacelle geometry is read as an STL triangulation (\c triSurface::New,
    auto-detecting ASCII/binary). One node per triangle is used: the node
    position is the triangle centroid, the node normal is the outward unit
    face normal, and the node area weight is the face area (design Decision 2).
    The node ordering is the surface face order, which is deterministic for a
    fixed input file.

    The force model is the paper's bluff-body closure with no blade-element
    (BEM) content:

    - normal force (Eq. 19, direct forcing, \c u^d = 0):
      \f$ f_n = h (\tilde{u} \cdot e_n) / \Delta t \, e_n \f$
    - tangential force (Eq. 21):
      \f$ f_\tau = \frac{1}{2} c_f U^2 e_\tau \f$
    - friction coefficient (Eq. 22, Schultz-Grunow) with an optional constant
      override
    - tangential direction (Eq. 23): the normalized velocity at the off-wall
      probe point \f$ X + h e_n \f$

    The unit system is the paper's: \c f_n and \c f_\tau are tractions per unit
    area **per unit density** (m²/s²). They become per-node forces only after
    multiplication by the node area (m⁴/s², \c nodeForces_); the SI contract
    exposed by \c forces() (newtons) is \c rhoRef_ times that list.

    Two per-node force lists are kept, sharing the same node ordering:

    - \c nodeForces_ — per-unit-density force \f$ (f_n + f_\tau)_i A_i \f$
      (m⁴/s²), consumed by the source's distribution step. Stored in the
      **global** mesh frame, because that is the frame the kernel spread acts
      in.
    - \c forces_ — SI force on the body (N), \c rhoRef_ * \c nodeForces_,
      reported in the **body** frame by \c forces() (the S3 FSI contract).

    The body frame is defined by the optional \c bodyOrigin and \c bodyAxis
    entries (defaults \c (0 0 0) and \c (1 0 0)); it defaults to the global
    mesh frame. \c bodyAxis is the nacelle axis expressed in the global frame,
    and the remaining body axes are completed deterministically. Positions,
    normals, and forces are transformed into the body frame by the accessors.

SourceFiles
    nacelleSurfaceSampler.C

\*---------------------------------------------------------------------------*/

#ifndef nacelleSurfaceSampler_H
#define nacelleSurfaceSampler_H

#include "List.H"
#include "dictionary.H"
#include "point.H"
#include "vector.H"
#include "tensor.H"
#include "triSurface.H"
#include "fvMesh.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

class nacelleSurfaceSource;

/*---------------------------------------------------------------------------*\
                    Class nacelleSurfaceSampler Declaration
\*---------------------------------------------------------------------------*/

class nacelleSurfaceSampler
{
    // Private data

        //- Allow the owning source to populate the per-node force lists
        friend class nacelleSurfaceSource;

        //- Mesh the surface is sampled on
        const fvMesh& mesh_;

        //- Surface triangulation, auto-detected ASCII/binary
        autoPtr<triSurface> surface_;

        //- Node positions (global frame, triangle centroids, face order)
        List<point> positions_;

        //- Node outward unit normals (global frame)
        List<vector> normals_;

        //- Node area weights
        List<scalar> areas_;

        //- Node positions in the body frame (the sampling contract)
        List<point> positionsBody_;

        //- Node normals in the body frame (the sampling contract)
        List<vector> normalsBody_;

        //- Per-node force ON the body, SI newtons, body frame (the contract)
        //  rhoRef_*(f_n + f_tau)_i*A_i
        List<vector> forces_;

        //- Per-node force on the body per unit density, global frame
        //  (f_n + f_tau)_i*A_i [m^4/s^2]; input to the distribution step
        List<vector> nodeForces_;

        //- Body origin in the global frame
        vector bodyOrigin_;

        //- Body-to-global rotation tensor (columns: body axes in global frame)
        tensor bodyToGlobal_;

        //- Body-to-global transform is identity
        bool identityBodyFrame_;

        //- Reference incoming velocity U (Eqs. 21, 22)
        scalar referenceVelocity_;

        //- Kinematic viscosity used by the Schultz-Grunow relation
        scalar nu_;

        //- Constant cf override; negative means use Schultz-Grunow
        scalar cfOverride_;

        //- Reference density for the SI-newton conversion of forces()
        scalar rhoRef_;

        //- Streamwise direction (the body axis, the nacelle axis is streamwise)
        vector streamwiseDirection_;

        //- Streamwise coordinate of the most-upstream node (the nose)
        scalar noseStreamwise_;

        //- Streamwise coordinate of each node
        List<scalar> streamwiseCoords_;


    // Private Member Functions

        //- No copy construct
        nacelleSurfaceSampler(const nacelleSurfaceSampler&) = delete;

        //- No copy assignment
        void operator=(const nacelleSurfaceSampler&) = delete;

        //- Build the body frame transform and the body-frame node geometry
        void createBodyFrame(const dictionary& dict);

        //- Resolve nu from the dictionary or constant/transportProperties
        void readViscosity(const dictionary& dict);


public:

    // Constructors

        //- Construct from a dictionary and the mesh
        nacelleSurfaceSampler(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    ~nacelleSurfaceSampler();


    // Member Functions

        // Access

            //- Number of nodes (triangle centroids)
            inline label nNodes() const;

            //- Node positions in the body frame
            inline const List<point>& positions() const;

            //- Per-node force on the body, SI newtons, body frame
            inline const List<vector>& forces() const;

            //- Node outward unit normals in the body frame
            inline const List<vector>& normals() const;

            //- Node area weights
            inline const List<scalar>& areas() const;

            //- Reference density used for the SI conversion
            inline scalar rhoRef() const;

            //- Reference incoming velocity U (Eqs. 21, 22)
            inline scalar referenceVelocity() const;

            //- Node positions in the global (mesh) frame
            inline const List<point>& positionsGlobal() const;

            //- Node outward unit normals in the global (mesh) frame
            inline const List<vector>& normalsGlobal() const;

            //- Body-to-global rotation tensor (columns: body axes in the
            //  global frame); its transpose maps global vectors into the
            //  body frame
            inline const tensor& bodyToGlobal() const;


        // Geometry

            //- Local cell size h = cbrt(V[cell]) of the cell containing X,
            //  resolved with a findCell + reduce(minOp) sentinel. Fatal if the
            //  point is not found on any processor (unreachable sample).
            scalar cellSize(const point& X) const;


        // Force model (per-node, per unit density)

            //- Normal traction (Eq. 19) with u^d = 0
            vector normalForce
            (
                const label i,
                const vector& uTilde,
                const scalar h,
                const scalar dt
            ) const;

            //- Tangential traction (Eq. 21) along eTau
            vector tangentialForce
            (
                const label i,
                const scalar cf,
                const vector& eTau
            ) const;

            //- Friction coefficient (Eq. 22), or the constant override
            scalar frictionCoefficient(const label i) const;

            //- Tangential direction (Eq. 23) at the off-wall probe point
            vector tangentialDirection
            (
                const label i,
                const volVectorField& U,
                const scalar h
            ) const;


        // Kernel and interpolation

            //- Smoothed four-point cosine kernel (Eq. 8), support |r| <= 2.5
            static scalar kernel(const scalar r);

            //- Kernel-smoothed velocity interpolation at X (Eq. 7); h is the
            //  node's local cell size
            vector interpolateVelocity
            (
                const point& X,
                const volVectorField& U,
                const scalar h
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nacelleSurfaceSamplerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
