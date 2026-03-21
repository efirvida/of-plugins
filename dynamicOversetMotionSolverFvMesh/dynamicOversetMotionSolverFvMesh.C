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
    defineTypeNameAndDebug(dynamicOversetMotionSolverFvMesh, 0);

    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicOversetMotionSolverFvMesh,
        IOobject
    );

    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicOversetMotionSolverFvMesh,
        doInit
    );
}


Foam::dynamicOversetMotionSolverFvMesh::dynamicOversetMotionSolverFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicMotionSolverFvMesh(io, doInit),
    oversetFvMeshBase(static_cast<const fvMesh&>(*this))
{}


Foam::dynamicOversetMotionSolverFvMesh::~dynamicOversetMotionSolverFvMesh()
{}


bool Foam::dynamicOversetMotionSolverFvMesh::update()
{
    if (!dynamicMotionSolverFvMesh::update())
    {
        return false;
    }

    oversetFvMeshBase::update();

    return true;
}


bool Foam::dynamicOversetMotionSolverFvMesh::writeObject
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
