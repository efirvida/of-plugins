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
    Foam::fv::bladeSurfaceSampler

Description
    Rotating-frame surface sampler for the blade actuator surface model
    (Yang & Sotiropoulos, arXiv:1702.02108v4, Sec. 2.1 extended to an imported
    surface; design D4-D7).

    The sampler reads the blade triangulation in the canonical generation
    frame of the geometry pipeline (+z = span outward, +y = reference chord
    trailing -> leading, +x = y x z; design section 4.4) and places it in the
    rotor frame through the construction frame injected by the owning turbine
    source:

        X_global = surfaceOrigin
                 + z_gen*surfaceSpanDirection
                 + y_gen*surfaceChordDirection
                 + x_gen*(surfaceChordDirection ^ surfaceSpanDirection)

    Every triangle centroid becomes a node carrying the S1 sampling assets
    (area weight, outward normal, stable surface order). Each node is then
    associated once with an actuator line element:

    - station:        s(X) = (X - surfaceOrigin) & surfaceSpanDirection
    - chord fraction: c(X) = chordMount_i
                             - ((X - position_i) & unit(chordDirection_i))
                               /chordLength_i
      (0 at the leading edge, 1 at the trailing edge; the element frame of the
      assigned patch, design section 4.4)

    The node-to-element partition is the 1-D Voronoi partition along the
    station (patch boundaries at the element-station midpoints, open at the
    root/tip ends). The constructor asserts, fatally, that every node is
    assigned, that every element owns a non-empty patch and that the patch
    areas sum to the total triSurface area.

    The distribution uses a per-node candidate list built once on top of a
    per-rank uniform bin grid over the local cell centres (bin size
    max(min_i support_i, small)). Each query visits the bins overlapping the
    node's kernel support and applies the S1 axis-aligned prefilter, so a
    per-addSup update is O(sum_i |candidates_i|) instead of a full local-cell
    scan. The support is 2.5*h_i for the paper cosine kernel (Eq. 8) or
    eps_i*sqrt(ln(1000)) with eps_i = 2*cbrt(V_i)*meshFactor for the Gaussian
    ablation (D7).

    Body frame contract: the frame is built from the optional bodyOrigin and
    bodyAxis entries (defaults surfaceOrigin and surfaceSpanDirection), i.e.
    fixed at construction, so positions()/normals() are azimuth-dependent
    while the node ordering is fixed. Every rigid transform the owning blade
    forwards (rotate/translate/pitch) moves the global nodes and refreshes the
    body-frame copies. pitch() is the documented rigid approximation about the
    root element pitch axis through the root chord pitch-axis point (D5).

    MPI follows the S1 replicated-node/local-cell design: every rank holds the
    full node list and builds candidates over its own local cells only, so
    every local cell is written exactly once by its owner.

SourceFiles
    bladeSurfaceSampler.C

\*---------------------------------------------------------------------------*/

#ifndef bladeSurfaceSampler_H
#define bladeSurfaceSampler_H

#include "surfaceSamplerBase.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

class actuatorLineElement;
class bladeSurfaceSource;

/*---------------------------------------------------------------------------*\
                  Class bladeSurfaceSampler Declaration
\*---------------------------------------------------------------------------*/

class bladeSurfaceSampler
:
    public surfaceSamplerBase
{
    // Private data

        //- Allow the owning distributor to read the partition and to fill the
        //  per-node force lists (the S1 friend pattern)
        friend class bladeSurfaceSource;

        //- Node positions in the canonical generation frame (raw STL order)
        List<point> canonicalNodes_;

        //- Span station of every node, m
        List<scalar> station_;

        //- LE-based chord fraction of every node (0 at the LE, 1 at the TE)
        List<scalar> chordFraction_;

        //- Owning element index of every node
        List<label> patch_;

        //- Node share of its patch force, A_i/sum_patch(A)
        List<scalar> patchAreaShare_;

        //- Local cell size of every node, cbrt(V) of its containing cell
        List<scalar> nodeH_;

        //- Per-node force on the blade per unit density, global frame
        //  patchAreaShare_[i]*elements_[patch_[i]].force()
        List<vector> nodeForces_;

        //- Per-node force on the blade in SI newtons, body frame (contract)
        List<vector> forces_;

        //- Bounded candidate local cell labels per node
        List<List<label>> candidates_;

        //- Span station of every element, m
        List<scalar> elementStation_;

        //- Total patch area of every element, m^2
        List<scalar> patchArea_;

        //- Rotor origin of the construction frame (global frame)
        vector surfaceOrigin_;

        //- Construction span direction, outward root -> tip (unit)
        vector spanAxis_;

        //- Construction chord direction, trailing -> leading (unit)
        vector chordAxis_;

        //- Construction section normal, chordAxis_ ^ spanAxis_ (unit)
        vector normalAxis_;

        //- Root element pitch axis (element spanDirection at construction)
        vector pitchAxis_;

        //- Root chord pitch-axis point (root quarter-chord at construction)
        point pitchPoint_;

        //- Kernel mode: cosine (paper Eq. 8) or gaussian (D7 ablation)
        word kernelMode_;

        //- Gaussian width factor (meshFactor), D7 fallback chain
        scalar meshFactor_;


    // Private Member Functions

        //- No copy construct
        bladeSurfaceSampler(const bladeSurfaceSampler&) = delete;

        //- No copy assignment
        void operator=(const bladeSurfaceSampler&) = delete;

        //- Place the canonical STL nodes through the construction frame
        void placeCanonicalNodes();

        //- Station/chord association and the 1-D Voronoi patch partition
        void buildPartition(PtrList<actuatorLineElement>& elements);

        //- Kernel support radius of a node
        scalar nodeSupport(const scalar h) const;

        //- Build the per-rank bin grid and the per-node candidate lists
        void buildCandidates();

        //- Refresh the body-frame positions/normals from the global lists
        void updateBodyFrame();

        //- Rotate the node geometry rigidly about point/axis and refresh the
        //  body frame
        void rotateGeometry
        (
            const point& rotationPoint,
            const vector& axis,
            const scalar radians
        );

        //- Rotation matrix about axis by radians; the same arithmetic as
        //  actuatorLineElement::rotate (from SOWFA), for lockstep
        static tensor rotationMatrix(const vector& axis, const scalar radians);

        //- Distribute the node forces onto their candidate cells (Eq. 18);
        //  the kernel is the paper cosine or the Gaussian ablation
        void distributeNodeForces
        (
            volVectorField& ff,
            const List<vector>& nodeForces
        ) const;


public:

    // Constructors

        //- Construct from the blade subdictionary, the mesh and the owning
        //  blade elements (the element frame drives the station/chord
        //  association and the pitch axis)
        bladeSurfaceSampler
        (
            const dictionary& dict,
            const fvMesh& mesh,
            PtrList<actuatorLineElement>& elements
        );


    //- Destructor
    virtual ~bladeSurfaceSampler();


    // Member Functions

        // Access

            //- Node positions in the canonical generation frame
            inline const List<point>& canonicalNodes() const;

            //- Span station of every node, m
            inline const List<scalar>& station() const;

            //- LE-based chord fraction of every node
            inline const List<scalar>& chordFraction() const;

            //- Owning element index of every node
            inline const List<label>& patch() const;

            //- Node share of its patch force
            inline const List<scalar>& patchAreaShare() const;

            //- Local cell size of every node
            inline const List<scalar>& nodeH() const;

            //- Span station of every element, m
            inline const List<scalar>& elementStation() const;

            //- Total patch area of every element, m^2
            inline const List<scalar>& patchArea() const;

            //- Per-node force on the blade in SI newtons, body frame
            inline const List<vector>& forces() const;

            //- Kernel mode ('cosine' or 'gaussian')
            inline const word& kernelMode() const;

            //- Gaussian width factor
            inline scalar meshFactor() const;


        // Manipulation (rigid, forwarded by the owning blade source)

            //- Rotate the nodes about point/axis by radians
            void rotate
            (
                const point& rotationPoint,
                const vector& axis,
                const scalar radians
            );

            //- Translate the nodes
            void translate(const vector& translation);

            //- Pitch the surface rigidly about the root element pitch axis
            //  through the root chord pitch-axis point (documented
            //  approximation, D5)
            void pitch(const scalar radians);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bladeSurfaceSamplerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
