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

#include "dynamicOversetZoneDisplacementFvMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(dynamicOversetZoneDisplacementFvMesh, 0);

    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicOversetZoneDisplacementFvMesh,
        IOobject
    );
}


Foam::dynamicOversetZoneDisplacementFvMesh::dynamicOversetZoneDisplacementFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    // Defer init: the motion solver created during init() may trigger
    // overset stencil construction (via inverseDistanceDiffusivity →
    // wallDist → oversetFvPatchField::initEvaluate).  The stencil
    // requires the full dynamicOversetZoneDisplacementFvMesh vtable
    // and oversetFvMeshBase to be constructed first.
    dynamicMotionSolverFvMesh(io, false),
    oversetFvMeshBase(static_cast<const fvMesh&>(*this))
{
    if (doInit)
    {
        if (!init(true))
        {
            FatalErrorInFunction
                << "Failed to initialise dynamicOversetZoneDisplacementFvMesh"
                << exit(FatalError);
        }
    }
}


bool Foam::dynamicOversetZoneDisplacementFvMesh::init(const bool doInit)
{
    // Force initial stencil computation.  The oversetFvMeshBase constructor
    // creates the cellCellStencilObject with update=false, leaving the
    // stencil empty (null cellInterpolationMap, empty cellTypes, etc.).
    // In the upstream dynamicOversetFvMesh the stencil is computed earlier
    // (during eagerly-created motion solvers) so oversetFvMeshBase finds
    // an already-updated stencil in the MeshObject cache.  Our deferred
    // init means oversetFvMeshBase creates the empty stencil first.
    // Without this update, any subsequent field evaluation (e.g. turbulence
    // model construction evaluating k/epsilon with overset BCs) triggers
    // oversetFvPatchField::initEvaluate() which accesses the null
    // cellInterpolationMap, corrupting the heap.
    //
    // The const_cast is safe here because:
    // - Stencil::New() returns a reference to the MeshObject-cached stencil
    // - movePoints() is a legitimate state update for the stencil
    // - The stencil is mutable by design; the const qualifier on New()
    //   reflects OpenFOAM's MeshObject interface, not immutability intent
    cellCellStencilObject& stencil =
        const_cast<cellCellStencilObject&>(Stencil::New(*this));

    stencil.movePoints();

    return dynamicMotionSolverFvMesh::init(doInit);
}


Foam::dynamicOversetZoneDisplacementFvMesh::~dynamicOversetZoneDisplacementFvMesh()
{
    // No explicit cleanup needed:
    // - oversetFvMeshBase destructor handles stencil via MeshObject lifecycle
    // - dynamicMotionSolverFvMesh destructor handles motion solver via autoPtr
    // - lduPtr_ is automatically cleaned up by autoPtr
}


bool Foam::dynamicOversetZoneDisplacementFvMesh::update()
{
    // Order is critical: motion must be applied before overset stencil update.
    // Reversing this order causes stale stencil data to be used during mesh
    // motion, leading to incorrect interpolation and potential solver failure.
    if (!dynamicMotionSolverFvMesh::update())
    {
        return false;
    }

    oversetFvMeshBase::update();

    return true;
}


bool Foam::dynamicOversetZoneDisplacementFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    bool ok1 =
        dynamicMotionSolverFvMesh::writeObject(streamOpt, writeOnProc);
    bool ok2 =
        oversetFvMeshBase::writeObject(streamOpt, writeOnProc);
    return ok1 && ok2;
}


// ************************************************************************* //
