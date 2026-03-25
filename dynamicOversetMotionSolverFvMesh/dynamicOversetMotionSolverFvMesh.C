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
    dynamicMotionSolverFvMesh(io, doInit),
    oversetFvMeshBase(static_cast<const fvMesh&>(*this))
{}


Foam::dynamicOversetZoneDisplacementFvMesh::~dynamicOversetZoneDisplacementFvMesh()
{
    // Release overset-owned addressing/interface caches before fvMesh teardown.
    oversetFvMeshBase::clearOut();
}


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
