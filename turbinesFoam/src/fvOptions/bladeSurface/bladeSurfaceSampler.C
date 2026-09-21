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

#include "bladeSurfaceSampler.H"
#include "actuatorLineElement.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

//- Base-sampler dictionary: the injected construction frame provides the
//  body-frame defaults; a user bodyOrigin/bodyAxis still wins (D5)
static dictionary bladeBaseDict(const dictionary& dict)
{
    dictionary baseDict(dict);

    if (!baseDict.found("bodyOrigin") && dict.found("surfaceOrigin"))
    {
        baseDict.add("bodyOrigin", dict.get<vector>("surfaceOrigin"));
    }

    if (!baseDict.found("bodyAxis") && dict.found("surfaceSpanDirection"))
    {
        baseDict.add("bodyAxis", dict.get<vector>("surfaceSpanDirection"));
    }

    return baseDict;
}

} // End namespace fv
} // End namespace Foam


void Foam::fv::bladeSurfaceSampler::placeCanonicalNodes()
{
    // Raw STL centroids in the canonical generation frame (design 4.4)
    canonicalNodes_ = positions_;

    // Columns are the construction axes expressed in the global frame:
    // X_global = surfaceOrigin + x_gen*normal + y_gen*chord + z_gen*span
    const tensor canonicalToGlobal
    (
        normalAxis_.x(), chordAxis_.x(), spanAxis_.x(),
        normalAxis_.y(), chordAxis_.y(), spanAxis_.y(),
        normalAxis_.z(), chordAxis_.z(), spanAxis_.z()
    );

    forAll(positions_, i)
    {
        positions_[i] = surfaceOrigin_ + (canonicalToGlobal & canonicalNodes_[i]);
    }

    forAll(normals_, i)
    {
        normals_[i] = canonicalToGlobal & normals_[i];
    }

    updateBodyFrame();
}


void Foam::fv::bladeSurfaceSampler::buildPartition
(
    PtrList<actuatorLineElement>& elements
)
{
    const label nElements = elements.size();

    // Element stations: the same monotone span coordinate as the nodes (D4)
    elementStation_.setSize(nElements);

    forAll(elements, e)
    {
        elementStation_[e] =
            (elements[e].position() - surfaceOrigin_) & spanAxis_;
    }

    station_.setSize(nNodes());

    forAll(positions_, i)
    {
        station_[i] = (positions_[i] - surfaceOrigin_) & spanAxis_;
    }

    // 1-D Voronoi partition: every node belongs to the element with the
    // nearest station, i.e. patch boundaries at the element-station midpoints
    // and open first/last patches at the root/tip ends
    patch_.setSize(nNodes());
    patchArea_.setSize(nElements, 0.0);

    forAll(station_, i)
    {
        scalar bestDist = VGREAT;
        label best = 0;

        forAll(elementStation_, e)
        {
            const scalar dist = mag(station_[i] - elementStation_[e]);

            if (dist < bestDist)
            {
                bestDist = dist;
                best = e;
            }
        }

        patch_[i] = best;
        patchArea_[best] += areas_[i];
    }

    // Fatal invariant: every element owns a non-empty patch (checked before
    // the area shares divide by the patch area)
    forAll(patchArea_, e)
    {
        if (patchArea_[e] <= 0.0)
        {
            FatalErrorInFunction
                << "Blade surface element patch " << e << " is empty: the "
                << "element station " << elementStation_[e]
                << " owns no triangle centroid" << nl
                << exit(FatalError);
        }
    }

    // Chord fraction from the assigned element's frame (LE-based, D4.4) and
    // the area share of each node in its patch (D4.3)
    chordFraction_.setSize(nNodes());
    patchAreaShare_.setSize(nNodes());

    forAll(positions_, i)
    {
        const label e = patch_[i];

        const vector chordUnit =
            elements[e].chordDirection()/mag(elements[e].chordDirection());

        chordFraction_[i] =
            elements[e].chordMount()
          - ((positions_[i] - elements[e].position()) & chordUnit)
            /elements[e].chordLength();

        patchAreaShare_[i] = areas_[i]/patchArea_[e];
    }

    // Fatal invariant: the patch areas sum to the total triSurface area
    scalar totalArea = 0.0;
    scalar patchAreaSum = 0.0;

    forAll(areas_, i)
    {
        totalArea += areas_[i];
    }

    forAll(patchArea_, e)
    {
        patchAreaSum += patchArea_[e];
    }

    if (mag(patchAreaSum - totalArea) > 1e-10*max(totalArea, VSMALL))
    {
        FatalErrorInFunction
            << "Blade surface patch areas (" << patchAreaSum
            << ") do not sum to the total surface area (" << totalArea << ")"
            << nl << exit(FatalError);
    }
}


Foam::scalar Foam::fv::bladeSurfaceSampler::nodeSupport(const scalar h) const
{
    if (kernelMode_ == "gaussian")
    {
        // Truncated Gaussian: eps = 2*cbrt(V)*meshFactor, support
        // eps*sqrt(ln(1000)) (actuatorSurfaceElement::applyForceField)
        return 2.0*h*meshFactor_*Foam::sqrt(Foam::log(1.0/0.001));
    }

    // Paper cosine kernel, support |r| <= 2.5 (Eq. 8)
    return 2.5*h;
}


void Foam::fv::bladeSurfaceSampler::buildCandidates()
{
    candidates_.setSize(nNodes());
    nodeH_.setSize(nNodes());

    scalar minSupport = VGREAT;

    forAll(positions_, i)
    {
        // Local cell size h = cbrt(V[cell]); fatal on an unreachable sample
        nodeH_[i] = cellSize(positions_[i]);
        minSupport = min(minSupport, nodeSupport(nodeH_[i]));
    }

    // A rank with no local cells owns no receiving cell: empty candidates
    if (mesh_.nCells() == 0)
    {
        return;
    }

    // Per-rank bounding box of the local cell centres
    point minC(VGREAT, VGREAT, VGREAT);
    point maxC(-VGREAT, -VGREAT, -VGREAT);

    forAll(mesh_.C(), cellI)
    {
        const point& C = mesh_.C()[cellI];

        minC.x() = min(minC.x(), C.x());
        minC.y() = min(minC.y(), C.y());
        minC.z() = min(minC.z(), C.z());
        maxC.x() = max(maxC.x(), C.x());
        maxC.y() = max(maxC.y(), C.y());
        maxC.z() = max(maxC.z(), C.z());
    }

    // Uniform bin grid, bin size max(min_i support_i, small) (D6)
    scalar binSize = max(minSupport, SMALL);

    label nBinsX = 0;
    label nBinsY = 0;
    label nBinsZ = 0;

    // Memory guard: a pathological support-to-box ratio must not allocate an
    // unbounded number of (mostly empty) bins. Enlarging the bins keeps the
    // query correct because the support prefilter below still removes every
    // cell outside the node's support box; it only adds a few candidate
    // entries per node.
    const scalar maxBins = 1e6;

    do
    {
        nBinsX = max(label((maxC.x() - minC.x())/binSize) + 1, label(1));
        nBinsY = max(label((maxC.y() - minC.y())/binSize) + 1, label(1));
        nBinsZ = max(label((maxC.z() - minC.z())/binSize) + 1, label(1));

        if (scalar(nBinsX)*scalar(nBinsY)*scalar(nBinsZ) <= maxBins)
        {
            break;
        }

        binSize *= 2.0;
    }
    while (true);

    List<List<label>> bins(nBinsX*nBinsY*nBinsZ);

    forAll(mesh_.C(), cellI)
    {
        const point& C = mesh_.C()[cellI];

        const label ix =
            min(max(label((C.x() - minC.x())/binSize), label(0)), nBinsX - 1);
        const label iy =
            min(max(label((C.y() - minC.y())/binSize), label(0)), nBinsY - 1);
        const label iz =
            min(max(label((C.z() - minC.z())/binSize), label(0)), nBinsZ - 1);

        bins[(iz*nBinsY + iy)*nBinsX + ix].append(cellI);
    }

    // Per-node query over the bins overlapping [X - support, X + support]
    forAll(positions_, i)
    {
        const point& X = positions_[i];
        const scalar s = nodeSupport(nodeH_[i]);

        const label i0x =
            min(max(label((X.x() - s - minC.x())/binSize), label(0)), nBinsX - 1);
        const label i1x =
            min(max(label((X.x() + s - minC.x())/binSize), label(0)), nBinsX - 1);
        const label i0y =
            min(max(label((X.y() - s - minC.y())/binSize), label(0)), nBinsY - 1);
        const label i1y =
            min(max(label((X.y() + s - minC.y())/binSize), label(0)), nBinsY - 1);
        const label i0z =
            min(max(label((X.z() - s - minC.z())/binSize), label(0)), nBinsZ - 1);
        const label i1z =
            min(max(label((X.z() + s - minC.z())/binSize), label(0)), nBinsZ - 1);

        for (label iz = i0z; iz <= i1z; ++iz)
        {
            for (label iy = i0y; iy <= i1y; ++iy)
            {
                for (label ix = i0x; ix <= i1x; ++ix)
                {
                    const List<label>& bin =
                        bins[(iz*nBinsY + iy)*nBinsX + ix];

                    forAll(bin, k)
                    {
                        const label cellI = bin[k];
                        const vector d = mesh_.C()[cellI] - X;

                        // S1 axis-aligned prefilter on the kernel support
                        if
                        (
                            mag(d.x()) > s
                         || mag(d.y()) > s
                         || mag(d.z()) > s
                        )
                        {
                            continue;
                        }

                        candidates_[i].append(cellI);
                    }
                }
            }
        }
    }
}


void Foam::fv::bladeSurfaceSampler::updateBodyFrame()
{
    if (identityBodyFrame_)
    {
        positionsBody_ = positions_;
        normalsBody_ = normals_;
        return;
    }

    const tensor globalToBody = bodyToGlobal_.T();

    forAll(positions_, i)
    {
        positionsBody_[i] = globalToBody & (positions_[i] - bodyOrigin_);
        normalsBody_[i] = globalToBody & normals_[i];
    }
}


Foam::tensor Foam::fv::bladeSurfaceSampler::rotationMatrix
(
    const vector& axis,
    const scalar radians
)
{
    // Same arithmetic as actuatorLineElement::rotate (from SOWFA) so the
    // surface stays in lockstep with the blade elements
    tensor RM;
    const scalar angle = radians;

    RM.xx() = sqr(axis.x()) + (1.0 - sqr(axis.x()))*cos(angle);
    RM.xy() = axis.x()*axis.y()*(1.0 - cos(angle)) - axis.z()*sin(angle);
    RM.xz() = axis.x()*axis.z()*(1.0 - cos(angle)) + axis.y()*sin(angle);
    RM.yx() = axis.x()*axis.y()*(1.0 - cos(angle)) + axis.z()*sin(angle);
    RM.yy() = sqr(axis.y()) + (1.0 - sqr(axis.y()))*cos(angle);
    RM.yz() = axis.y()*axis.z()*(1.0 - cos(angle)) - axis.x()*sin(angle);
    RM.zx() = axis.x()*axis.z()*(1.0 - cos(angle)) - axis.y()*sin(angle);
    RM.zy() = axis.y()*axis.z()*(1.0 - cos(angle)) + axis.x()*sin(angle);
    RM.zz() = sqr(axis.z()) + (1.0 - sqr(axis.z()))*cos(angle);

    return RM;
}


void Foam::fv::bladeSurfaceSampler::rotateGeometry
(
    const point& rotationPoint,
    const vector& axis,
    const scalar radians
)
{
    const tensor RM = rotationMatrix(axis, radians);

    forAll(positions_, i)
    {
        positions_[i] = rotationPoint + (RM & (positions_[i] - rotationPoint));
    }

    forAll(normals_, i)
    {
        normals_[i] = RM & normals_[i];
    }

    // The pitch frame is rigidly attached to the surface
    pitchPoint_ = rotationPoint + (RM & (pitchPoint_ - rotationPoint));
    pitchAxis_ = RM & pitchAxis_;

    updateBodyFrame();
}


void Foam::fv::bladeSurfaceSampler::distributeNodeForces
(
    volVectorField& ff,
    const List<vector>& nodeForces
) const
{
    const scalarField& V = mesh_.V();

    forAll(nodeForces, i)
    {
        const point& X = positions_[i];
        const vector F = nodeForces[i];
        const scalar h = nodeH_[i];
        const List<label>& cells = candidates_[i];

        const scalar eps = 2.0*h*meshFactor_;

        forAll(cells, k)
        {
            const label cellI = cells[k];
            const vector d = mesh_.C()[cellI] - X;

            scalar w = 0.0;

            if (kernelMode_ == "gaussian")
            {
                // No-mesh ASM Gaussian ablation: the dimensionless cell
                // weight is V_c*phi so that F*(w/V_c) below is the element's
                // exp(-d^2/eps^2)/(eps^3*pi^1.5) factor
                w =
                    V[cellI]*Foam::exp(-magSqr(d)/sqr(eps))
                   /(Foam::pow(eps, 3)
                     *Foam::pow(constant::mathematical::pi, 1.5));
            }
            else
            {
                w = kernel(d.x()/h)*kernel(d.y()/h)*kernel(d.z()/h);
            }

            if (w > 0.0)
            {
                // Eq. 18 with the paper's negative sign: the distributed
                // force acts on the flow, F acts on the blade
                ff[cellI] -= F*(w/V[cellI]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::bladeSurfaceSampler::bladeSurfaceSampler
(
    const dictionary& dict,
    const fvMesh& mesh,
    PtrList<actuatorLineElement>& elements
)
:
    surfaceSamplerBase(bladeBaseDict(dict), mesh, "surfaceGeometry"),
    canonicalNodes_(),
    station_(),
    chordFraction_(),
    patch_(),
    patchAreaShare_(),
    nodeH_(),
    nodeForces_(),
    forces_(),
    candidates_(),
    elementStation_(),
    patchArea_(),
    surfaceOrigin_(dict.lookupOrDefault("surfaceOrigin", bodyOrigin_)),
    spanAxis_(vector(0, 0, 1)),
    chordAxis_(vector(0, 1, 0)),
    normalAxis_(vector(1, 0, 0)),
    pitchAxis_(vector(0, 0, 1)),
    pitchPoint_(vector::zero),
    kernelMode_(dict.lookupOrDefault<word>("kernel", "cosine")),
    meshFactor_(1.0)
{
    if (elements.size() == 0)
    {
        FatalErrorInFunction
            << "The blade surface sampler requires at least one actuator line "
            << "element to associate the surface nodes with" << nl
            << exit(FatalError);
    }

    if (kernelMode_ != "cosine" && kernelMode_ != "gaussian")
    {
        FatalErrorInFunction
            << "Unknown blade surface kernel '" << kernelMode_ << "'" << nl
            << "Valid kernels: cosine (paper Eq. 8), gaussian (ablation)"
            << nl << exit(FatalError);
    }

    // Construction frame (D5): the injected surface* entries are the defaults,
    // a user bodyOrigin/bodyAxis override wins
    spanAxis_ = dict.lookupOrDefault
    (
        "surfaceSpanDirection",
        dict.lookupOrDefault("bodyAxis", vector(0, 0, 1))
    );

    if (mag(spanAxis_) < SMALL)
    {
        FatalErrorInFunction
            << "The construction span direction must be non-zero" << nl
            << exit(FatalError);
    }
    spanAxis_ /= mag(spanAxis_);

    chordAxis_ = dict.lookupOrDefault
    (
        "surfaceChordDirection",
        vector(0, 1, 0)
    );

    if (mag(chordAxis_) < SMALL)
    {
        FatalErrorInFunction
            << "The construction chord direction must be non-zero" << nl
            << exit(FatalError);
    }
    chordAxis_ /= mag(chordAxis_);

    normalAxis_ = chordAxis_ ^ spanAxis_;

    if (mag(normalAxis_) < SMALL)
    {
        FatalErrorInFunction
            << "The construction span and chord directions must not be "
            << "parallel" << nl << exit(FatalError);
    }
    normalAxis_ /= mag(normalAxis_);

    // Pitch frame: the root element's own pitch axis and quarter-chord point
    // (D5; the element pitch axis is anti-parallel to the outward span in the
    // Phase VI orientation, so the raw radians are forwarded about the element
    // axis, not the construction span)
    pitchAxis_ = elements[0].spanDirection();

    if (mag(pitchAxis_) < SMALL)
    {
        FatalErrorInFunction
            << "The root element span direction must be non-zero" << nl
            << exit(FatalError);
    }
    pitchAxis_ /= mag(pitchAxis_);

    {
        const vector chordUnit =
            elements[0].chordDirection()/mag(elements[0].chordDirection());

        pitchPoint_ =
            elements[0].position()
          + (elements[0].chordMount() - 0.25)
           *elements[0].chordLength()*chordUnit;
    }

    // Gaussian width: the surface subdict, then the blade's profileData
    // GaussianCoeffs (the owning element dicts carry it), else the element
    // default 2.0 (D7)
    if (kernelMode_ == "gaussian")
    {
        meshFactor_ = dict.lookupOrDefault<scalar>("meshFactor", -1.0);

        if (meshFactor_ < 0.0)
        {
            const dictionary gauss =
                dict.subOrEmptyDict("profileData").subOrEmptyDict
                (
                    "GaussianCoeffs"
                );

            meshFactor_ = gauss.lookupOrDefault<scalar>("meshFactor", -1.0);
        }

        if (meshFactor_ < 0.0)
        {
            const dictionary gauss =
                elements[0].profileDict().subOrEmptyDict("GaussianCoeffs");

            meshFactor_ = gauss.lookupOrDefault<scalar>("meshFactor", 2.0);
        }

        if (meshFactor_ <= 0.0)
        {
            FatalErrorInFunction
                << "The Gaussian 'meshFactor' must be positive" << nl
                << exit(FatalError);
        }
    }

    // Canonical STL nodes -> rotor frame, station/chord association and the
    // patch partition
    placeCanonicalNodes();
    buildPartition(elements);

    nodeForces_.setSize(nNodes(), vector::zero);
    forces_.setSize(nNodes(), vector::zero);

    // Per-node cell size, kernel support and bounded candidate lists (once:
    // the background mesh and h_i are static)
    buildCandidates();

    Info<< "Blade surface sampler: " << nNodes() << " nodes over "
        << elements.size() << " element patches, kernel " << kernelMode_;

    if (kernelMode_ == "gaussian")
    {
        Info<< " (meshFactor " << meshFactor_ << ")";
    }

    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::bladeSurfaceSampler::~bladeSurfaceSampler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::bladeSurfaceSampler::rotate
(
    const point& rotationPoint,
    const vector& axis,
    const scalar radians
)
{
    if (mag(axis) < SMALL)
    {
        FatalErrorInFunction
            << "The rotation axis must be non-zero" << nl
            << exit(FatalError);
    }

    rotateGeometry(rotationPoint, axis/mag(axis), radians);
}


void Foam::fv::bladeSurfaceSampler::translate(const vector& translation)
{
    forAll(positions_, i)
    {
        positions_[i] += translation;
    }

    pitchPoint_ += translation;

    updateBodyFrame();
}


void Foam::fv::bladeSurfaceSampler::pitch(const scalar radians)
{
    // Documented rigid approximation (D5): the surface pitches about the root
    // element pitch axis through the root chord pitch-axis point
    rotateGeometry(pitchPoint_, pitchAxis_, radians);
}


// ************************************************************************* //
