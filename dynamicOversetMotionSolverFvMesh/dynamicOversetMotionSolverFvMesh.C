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

#include "dynamicOversetMotionSolverFvMesh.H"
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
        init(true);
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
    const_cast<cellCellStencilObject&>
    (
        Stencil::New(*this)
    ).movePoints();

    return dynamicMotionSolverFvMesh::init(doInit);
}


Foam::dynamicOversetZoneDisplacementFvMesh::~dynamicOversetZoneDisplacementFvMesh()
{}


bool Foam::dynamicOversetZoneDisplacementFvMesh::update()
{
    if (!dynamicMotionSolverFvMesh::update())
    {
        return false;
    }

    // Keep same ordering used by upstream dynamicOversetFvMesh:
    // first move mesh, then refresh overset addressing/stencil state.
    oversetFvMeshBase::update();

    return true;
}


bool Foam::dynamicOversetZoneDisplacementFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    bool ok =
        dynamicMotionSolverFvMesh::writeObject(streamOpt, writeOnProc);
    ok = oversetFvMeshBase::writeObject(streamOpt, writeOnProc) && ok;
    return ok;
}


// ************************************************************************* //
