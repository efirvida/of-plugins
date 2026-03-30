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
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "mapPolyMesh.H"
#include "solidBodyMotionFunction.H"
#include "transformField.H"
#include "fvOptions.H"
#include "cellSet.H"
#include "syncTools.H"
#include "DynamicList.H"

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


Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::
solidBodyDisplacementLaplacianZoneFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    SBMFPtr_(solidBodyMotionFunction::New(coeffDict(), mesh.time())),
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
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    pointLocation_(nullptr),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
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

        DynamicList<label> ptIDs(fvMesh_.nPoints());
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

    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< type() << ":" << nl
            << "    diffusivity       : " << diffusivityPtr_().type() << nl
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
    SBMFPtr_(solidBodyMotionFunction::New(coeffDict(), mesh.time())),
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
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    pointLocation_(nullptr),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
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
        const cellSet set(fvMesh_, cellSetName);
        cellIDs = set.toc();
    }

    const label nCells = returnReduce(cellIDs.size(), sumOp<label>());
    moveAllCells_ = (nCells == 0);

    if (!moveAllCells_)
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

        DynamicList<label> ptIDs(fvMesh_.nPoints());
        forAll(movePts, pointi)
        {
            if (movePts[pointi])
            {
                ptIDs.append(pointi);
            }
        }

        pointIDs_.transfer(ptIDs);
    }

    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::MUST_READ,
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
{}


Foam::motionDiffusivity&
Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::diffusivity()
{
    if (!diffusivityPtr_)
    {
        diffusivityPtr_ = motionDiffusivity::New
        (
            fvMesh_,
            coeffDict().lookup("diffusivity")
        );
    }

    return *diffusivityPtr_;
}


Foam::tmp<Foam::pointField>
Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::curPoints() const
{
    // FSI-DISP-LOG: log cellDisplacement and pointDisplacement BEFORE interpolate
    {
        const scalarField magCellInt(mag(cellDisplacement_.primitiveField()));
        const scalarField magPointPrim(mag(pointDisplacement_.primitiveField()));
        Info << "[FSI-DISP-LOG] curPoints BEFORE interpolate:"
             << " max|cellDisp.internal|=" << gMax(magCellInt)
             << " avg|cellDisp.internal|=" << gAverage(magCellInt)
             << " max|pointDisp.prim|=" << gMax(magPointPrim)
             << endl;
    }

    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    // FSI-DISP-LOG: log pointDisplacement after cell-to-point interpolation
    {
        const scalarField magPrim(mag(pointDisplacement_.primitiveField()));
        Info << "[FSI-DISP-LOG] curPoints after interpolate:"
             << " max|pointDisp.prim|=" << gMax(magPrim)
             << " avg|pointDisp.prim|=" << gAverage(magPrim)
             << endl;
    }

    if
    (
        pointDisplacement_.boundaryField().size()
     != cellDisplacement_.boundaryField().size()
    )
    {
        pointDisplacement_.correctBoundaryConditions();
    }

    tmp<pointField> tnewPoints(new pointField(points0()));

    if (moveAllCells_)
    {
        tnewPoints = transformPoints(SBMFPtr_().transformation(), points0());
    }
    else
    {
        pointField& transformedPts = tnewPoints.ref();

        UIndirectList<point>(transformedPts, pointIDs_) = transformPoints
        (
            SBMFPtr_().transformation(),
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

        twoDCorrectPoints(curPoints);

        // FSI-DISP-LOG: log max displacement from reference points0()
        {
            const pointField disp = curPoints - points0();
            Info << "[FSI-DISP-LOG] curPoints mesh motion:"
                 << " maxMag(curPoints - points0)=" << gMax(mag(disp))
                 << endl;
        }

        return tcurPoints;
    }
}


void Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::solve()
{
    // FSI-DISP-LOG: state at solve() entry (before anything)
    Info << "[FSI-DISP-LOG] solve() ENTRY max|cellDisp.internal|="
         << gMax(mag(cellDisplacement_.primitiveField())) << endl;

    movePoints(fvMesh_.points());

    diffusivity().correct();

    // FSI-DISP-LOG: after diffusivity().correct()
    Info << "[FSI-DISP-LOG] solve() after diffusivity().correct() max|cellDisp.internal|="
         << gMax(mag(cellDisplacement_.primitiveField()))
         << endl;

    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    fv::options& fvOptions(fv::options::New(fvMesh_));

    // Update boundary coefficients before building the laplacian matrix.
    // We use the same pattern as the standard displacementLaplacianFvMotionSolver.
    // On overset meshes, calling correctBoundaryConditions() would trigger
    // oversetFvPatchField::initEvaluate() which performs field interpolation
    // via mapDistribute, introducing unwanted MPI communication that can
    // corrupt heap metadata (detected as "corrupted double-linked list"
    // during the next inverseDistance::update() call).
    // On cyclicAMI meshes (without overset), correctBoundaryConditions()
    // is needed to prevent MPI deadlock from stale AMI data during the
    // non-orthogonal correction (fvc::grad inside fvm::laplacian).
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

    // Zero internal field before Laplacian solve to prevent hole-cell
    // singularity at time-window transitions. Overset hole cells have
    // zero diagonal in the Laplacian matrix; a non-zero initial guess
    // (from checkpoint restore) causes the smoother to divide by ~0,
    // producing overflow. The BC (set absolutely by Displacement::read)
    // drives the solution, so the initial guess does not affect correctness.
    cellDisplacement_.primitiveFieldRef() = vector::zero;

    // FSI-DISP-LOG: after updateCoeffs, before Laplacian matrix assembly
    Info << "[FSI-DISP-LOG] solve() before TEqn max|cellDisp.internal|="
         << gMax(mag(cellDisplacement_.primitiveField())) << endl;

    fvVectorMatrix TEqn
    (
        fvm::laplacian
        (
            dimensionedScalar("viscosity", dimViscosity, 1.0)
           *diffusivity().operator()(),
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )
     ==
        fvOptions(cellDisplacement_)
    );

    fvOptions.constrain(TEqn);
    TEqn.solveSegregatedOrCoupled();

    // FSI-DISP-LOG: after Laplacian solve, before fvOptions.correct()
    Info << "[FSI-DISP-LOG] solve() after TEqn.solve max|cellDisp.internal|="
         << gMax(mag(cellDisplacement_.primitiveField())) << endl;

    fvOptions.correct(cellDisplacement_);

    // FSI-DISP-LOG: after fvOptions.correct()
    Info << "[FSI-DISP-LOG] solve() after fvOptions.correct max|cellDisp.internal|="
         << gMax(mag(cellDisplacement_.primitiveField())) << endl;
}


void Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);

    diffusivityPtr_.clear();
}


// ************************************************************************* //
