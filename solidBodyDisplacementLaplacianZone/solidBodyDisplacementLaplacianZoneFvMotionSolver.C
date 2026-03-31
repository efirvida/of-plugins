/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM plugin extensions.

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

#include "solidBodyDisplacementLaplacianZoneFvMotionSolver.H"

#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "IStringStream.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "mapPolyMesh.H"
#include "solidBodyMotionFunction.H"
#include "transformField.H"
#include "fvOptions.H"
#include "cellSet.H"
#include "syncTools.H"
#include "DynamicList.H"
#include "pointMesh.H"

namespace Foam
{
    defineTypeNameAndDebug(solidBodyDisplacementLaplacianZoneFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        solidBodyDisplacementLaplacianZoneFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        solidBodyDisplacementLaplacianZoneFvMotionSolver,
        displacement
    );
}


void Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::initZoneAndPointIDs()
{
    const word cellZoneName(coeffDict().getOrDefault<word>("cellZone", "none"));
    const word cellSetName(coeffDict().getOrDefault<word>("cellSet", "none"));

    if ((cellZoneName != "none") && (cellSetName != "none"))
    {
        FatalIOErrorInFunction(coeffDict())
            << "Either cellZone OR cellSet can be supplied, but not both. "
            << "If neither is supplied, all cells will be included"
            << exit(FatalIOError);
    }

    labelList cellIDs;

    if (cellZoneName != "none")
    {
        Info<< "Applying solid body motion to cellZone " << cellZoneName
            << endl;

        const label zoneID = fvMesh_.cellZones().findZoneID(cellZoneName);

        if (zoneID == -1)
        {
            FatalErrorInFunction
                << "Unable to find cellZone " << cellZoneName
                << ". Valid cellZones are: "
                << fvMesh_.cellZones().names()
                << exit(FatalError);
        }

        cellIDs = fvMesh_.cellZones()[zoneID];
    }

    if (cellSetName != "none")
    {
        Info<< "Applying solid body motion to cellSet " << cellSetName
            << endl;

        const cellSet set(fvMesh_, cellSetName);
        cellIDs = set.toc();
    }

    const label nCells = returnReduce(cellIDs.size(), sumOp<label>());
    moveAllCells_ = (nCells == 0);

    if (moveAllCells_)
    {
        Info<< "Applying rigid transform to entire mesh" << endl;
    }
    else
    {
        boolList movePts(fvMesh_.nPoints(), false);

        forAll(cellIDs, i)
        {
            const label celli = cellIDs[i];
            const cell& c = fvMesh_.cells()[celli];

            forAll(c, j)
            {
                const face& f = fvMesh_.faces()[c[j]];
                forAll(f, k)
                {
                    movePts[f[k]] = true;
                }
            }
        }

        syncTools::syncPointList(fvMesh_, movePts, orEqOp<bool>(), false);

        DynamicList<label> ptIDs;
        forAll(movePts, pointi)
        {
            if (movePts[pointi])
            {
                ptIDs.append(pointi);
            }
        }

        pointIDs_.transfer(ptIDs);

        Info<< "Applying rigid transform to "
            << returnReduce(pointIDs_.size(), sumOp<label>())
            << " points from selected region" << endl;
    }
}


Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::
solidBodyDisplacementLaplacianZoneFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    solidBodyMotionPtr_(solidBodyMotionFunction::New(coeffDict(), mesh.time())),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(dimLength, Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement().boundaryField())
    ),
    pointLocation_(nullptr),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, IStringStream(coeffDict().get<word>("interpolation"))())
      : motionInterpolation::New(fvMesh_)
    ),
    // Defer diffusivity creation to the first call to diffusivity().
    // Eager construction triggers wallDist → boundary evaluation which
    // can crash on overset meshes (oversetFvPatchField::initEvaluate
    // builds the cellCellStencil before the mesh object is fully set up).
    diffusivityPtr_(nullptr),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    ),
    pointIDs_(),
    moveAllCells_(false)
{
    if (!solidBodyMotionPtr_)
    {
        FatalErrorInFunction
            << "Failed to create solidBodyMotionFunction. "
            << "Check solidBodyMotionFunction entry in dictionary."
            << exit(FatalError);
    }

    initZoneAndPointIDs();

    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< type() << ":" << nl
            << "    diffusivity       : "
            << coeffDict().get<word>("diffusivity") << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
    }

    if (io.typeHeaderOk<pointVectorField>(true))
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMesh_)
            )
        );

        if (debug)
        {
            Info<< type() << " : Read pointVectorField " << io.name()
                << " for point boundary conditions. Boundary conditions: "
                << pointLocation_().boundaryField().types() << endl;
        }
    }
}


Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::
solidBodyDisplacementLaplacianZoneFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    solidBodyMotionPtr_(solidBodyMotionFunction::New(coeffDict(), mesh.time())),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(dimLength, Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement.boundaryField())
    ),
    pointLocation_(nullptr),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, IStringStream(coeffDict().get<word>("interpolation"))())
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_(nullptr),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    ),
    pointIDs_(),
    moveAllCells_(false)
{
    if (!solidBodyMotionPtr_)
    {
        FatalErrorInFunction
            << "Failed to create solidBodyMotionFunction. "
            << "Check solidBodyMotionFunction entry in dictionary."
            << exit(FatalError);
    }

    initZoneAndPointIDs();

    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (io.typeHeaderOk<pointVectorField>(true))
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMesh_)
            )
        );
    }
}


Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::
~solidBodyDisplacementLaplacianZoneFvMotionSolver()
{
    diffusivityPtr_.clear();
}


Foam::motionDiffusivity&
Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::diffusivity()
{
    if (!diffusivityPtr_)
    {
        diffusivityPtr_ = motionDiffusivity::New
        (
            fvMesh_,
            IStringStream(coeffDict().get<word>("diffusivity"))()
        );
    }

    return *diffusivityPtr_;
}


Foam::tmp<Foam::pointField>
Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    pointDisplacement_.correctBoundaryConditions();

    tmp<pointField> tnewPoints(new pointField(points0()));

    if (moveAllCells_)
    {
        tnewPoints = transformPoints
        (
            solidBodyMotionPtr_().transformation(),
            points0()
        );
    }
    else
    {
        pointField& transformedPts = tnewPoints.ref();

        UIndirectList<point>(transformedPts, pointIDs_) = transformPoints
        (
            solidBodyMotionPtr_().transformation(),
            pointField(points0(), pointIDs_)
        );
    }

    const pointField& newPoints = tnewPoints();

    if (pointLocation_)
    {
        pointLocation_().primitiveFieldRef() =
            newPoints + pointDisplacement_.internalField();

        pointLocation_().correctBoundaryConditions();

        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                pointLocation_()[pz[i]] = newPoints[pz[i]];
            }
        }

        // 2D correction: enforce zero displacement in third dimension
        // for 2D axisymmetric or wedge geometries
        twoDCorrectPoints(pointLocation_().primitiveFieldRef());

        return tmp<pointField>(pointLocation_().primitiveField());
    }
    else
    {
        tmp<pointField> tcurPoints
        (
            newPoints + pointDisplacement_.primitiveField()
        );
        pointField& curPoints = tcurPoints.ref();

        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                curPoints[pz[i]] = newPoints[pz[i]];
            }
        }

        // 2D correction: enforce zero displacement in third dimension
        // for 2D axisymmetric or wedge geometries
        twoDCorrectPoints(curPoints);

        // Reset cellDisplacement_ to zero after the mesh motion has been
        // applied.  This field is only a scratch variable used inside solve()
        // to drive the Laplacian; once curPoints() has consumed it, it has no
        // physical meaning. Zeroing here ensures that the preCICE checkpoint
        // (written immediately after curPoints() returns) always stores a clean
        // zero field.  Without this, the checkpoint carries the Laplacian
        // solution from the previous window (~8.6e-10), which becomes the
        // initial guess for the MOVED-MESH Laplacian in the next window and
        // causes overflow due to changed Laplacian coefficients.
        cellDisplacement_.primitiveFieldRef() = Zero;
        cellDisplacement_.correctBoundaryConditions();

        return tcurPoints;
    }
}


void Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::solve()
{
    movePoints(fvMesh_.points());

    diffusivity().correct();

    // Zero internal field BEFORE updateCoeffs so that the overset patch
    // interpolation (which reads donor cell values from the internal field)
    // starts from zero in every time-window iteration.  Without this, the
    // preCICE checkpoint restore leaves a non-zero internal field at the
    // start of window 2+, which the overset BC picks up as donor values,
    // creating a non-trivial initial condition that interacts badly with
    // the near-singular hole-cell rows in the Laplacian matrix.
    cellDisplacement_.primitiveFieldRef() = vector::zero;

    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    fv::options& fvOptions(fv::options::New(fvMesh_));

    // On overset meshes use updateCoeffs (not correctBoundaryConditions)
    // to avoid oversetFvPatchField::initEvaluate triggering mapDistribute
    // MPI comms that corrupt heap metadata.
    bool hasOversetPatch = false;
    forAll(cellDisplacement_.boundaryField(), patchi)
    {
        if (cellDisplacement_.boundaryField()[patchi].type() == "overset")
        {
            hasOversetPatch = true;
            break;
        }
    }

    if (hasOversetPatch)
    {
        cellDisplacement_.boundaryFieldRef().updateCoeffs();
    }
    else
    {
        cellDisplacement_.correctBoundaryConditions();
    }

    fvVectorMatrix TEqn
    (
        fvm::laplacian
        (
            dimensionedScalar("unitScale", dimless, 1.0)
           *diffusivity().operator()(),
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )
     ==
        fvOptions(cellDisplacement_)
    );

    // Constrain hole cells to zero displacement by adding a large diagonal
    // contribution.  Interior hole cells (cellMask < holeCellThreshold_) are
    // not on any patch so oversetFvPatchField::manipulateMatrix() does not
    // constrain them.  Their Laplacian diagonal can be near-zero (all
    // neighbours are also hole cells), causing the Gauss-Seidel smoother to
    // divide by ~0 and produce overflow.  Pinning them to zero is correct:
    // the motion field has no physical meaning inside the overset hole region.
    const volScalarField* cellMaskPtr =
        fvMesh_.findObject<volScalarField>("cellMask");
    if (cellMaskPtr)
    {
        const scalarField& mask = cellMaskPtr->primitiveField();
        scalarField& diag = TEqn.diag();
        forAll(mask, celli)
        {
            if (mask[celli] < holeCellThreshold_)
                diag[celli] += VGREAT;
        }
    }

    fvOptions.constrain(TEqn);
    TEqn.solveSegregatedOrCoupled();

    fvOptions.correct(cellDisplacement_);

    // Post-solve: zero hole cells using cellMask as a safety net.
    if (cellMaskPtr)
    {
        const scalarField& mask = cellMaskPtr->primitiveField();
        vectorField& dispRef = cellDisplacement_.primitiveFieldRef();
        forAll(dispRef, celli)
        {
            if (mask[celli] < holeCellThreshold_)
                dispRef[celli] = Zero;
        }
    }

}


void Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);

    diffusivityPtr_.clear();

    // Rebuild point IDs after topology change (refinement, layer addition)
    // since pointIDs_ may reference stale indices
    if (mpm.hasMotionPoints())
    {
        initZoneAndPointIDs();
    }
}


// ************************************************************************* //
