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

Class
    Foam::dynamicOversetZoneDisplacementFvMesh

Description
    Overset-enabled dynamicFvMesh wrapper for a single motion solver.

    This class follows the overset integration pattern from OpenFOAM's native
    dynamicOversetFvMesh, but uses dynamicMotionSolverFvMesh (single motion
    solver) instead of dynamicMotionSolverListFvMesh (multibody list).

    The goal is to keep robust overset lifecycle/update/write behavior while
    preserving the single-solver workflow used in FSI cases with rotating-zone
    motion plus displacement-based mesh deformation.

    Typical usage:

    \verbatim
    dynamicFvMesh   dynamicOversetZoneDisplacementFvMesh;

    motionSolverLibs
    (
        "libsolidBodyDisplacementLaplacianZoneFvMotionSolver.so"
        "libfvMotionSolvers.so"
        "libfsiOmega.so"
    );

    motionSolver    solidBodyDisplacementLaplacianZone;

    solidBodyDisplacementLaplacianZoneCoeffs
    {
        cellZone rotor;
        solidBodyMotionFunction rotatingMotion;
        rotatingMotionCoeffs
        {
            origin  (0 0 0);
            axis    (0 1 0);
            omega   158;
        }
        diffusivity inverseDistance (propellerTip);
    }
    \endverbatim

SourceFiles
    dynamicOversetMotionSolverFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicOversetMotionSolverFvMesh_H
#define dynamicOversetMotionSolverFvMesh_H

#include "motionSolver.H"
#include "dynamicMotionSolverFvMesh.H"
#include "oversetFvMeshBase.H"
#include "cellCellStencilObject.H"
#include "fvMatrices.H"
#include "volFields.H"

namespace Foam
{

class dynamicOversetZoneDisplacementFvMesh
:
    public dynamicMotionSolverFvMesh,
    public oversetFvMeshBase
{
    // Private Member Functions

        dynamicOversetZoneDisplacementFvMesh
        (
            const dynamicOversetZoneDisplacementFvMesh&
        ) = delete;

        void operator=(const dynamicOversetZoneDisplacementFvMesh&) = delete;


public:

    TypeName("dynamicOversetZoneDisplacementFvMesh");


    dynamicOversetZoneDisplacementFvMesh
    (
        const IOobject& io,
        const bool doInit=true
    );


    virtual ~dynamicOversetZoneDisplacementFvMesh();


    //- Deferred initialisation (called after full construction)
    virtual bool init(const bool doInit);


    virtual const lduAddressing& lduAddr() const
    {
        return oversetFvMeshBase::lduAddr();
    }

    virtual lduInterfacePtrsList interfaces() const
    {
        return oversetFvMeshBase::interfaces();
    }


    virtual void interpolate(scalarField& psi) const
    {
        Stencil::New(*this).interpolate<scalar>(*this, psi);
    }

    virtual void interpolate(vectorField& psi) const
    {
        Stencil::New(*this).interpolate<vector>(*this, psi);
    }

    virtual void interpolate(sphericalTensorField& psi) const
    {
        Stencil::New(*this).interpolate<sphericalTensor>(*this, psi);
    }

    virtual void interpolate(symmTensorField& psi) const
    {
        Stencil::New(*this).interpolate<symmTensor>(*this, psi);
    }

    virtual void interpolate(tensorField& psi) const
    {
        Stencil::New(*this).interpolate<tensor>(*this, psi);
    }

    virtual void interpolate(volScalarField& psi) const
    {
        Stencil::New(*this).interpolate<volScalarField>(psi);
    }

    virtual void interpolate(volVectorField& psi) const
    {
        Stencil::New(*this).interpolate<volVectorField>(psi);
    }

    virtual void interpolate(volSphericalTensorField& psi) const
    {
        Stencil::New(*this).interpolate<volSphericalTensorField>(psi);
    }

    virtual void interpolate(volSymmTensorField& psi) const
    {
        Stencil::New(*this).interpolate<volSymmTensorField>(psi);
    }

    virtual void interpolate(volTensorField& psi) const
    {
        Stencil::New(*this).interpolate<volTensorField>(psi);
    }


    virtual SolverPerformance<scalar> solve
    (
        fvMatrix<scalar>& m,
        const dictionary& dict
    ) const
    {
        return oversetFvMeshBase::solveOverset<scalar>(m, dict);
    }

    virtual SolverPerformance<vector> solve
    (
        fvMatrix<vector>& m,
        const dictionary& dict
    ) const
    {
        return oversetFvMeshBase::solveOverset<vector>(m, dict);
    }

    virtual SolverPerformance<symmTensor> solve
    (
        fvMatrix<symmTensor>& m,
        const dictionary& dict
    ) const
    {
        return oversetFvMeshBase::solveOverset<symmTensor>(m, dict);
    }

    virtual SolverPerformance<tensor> solve
    (
        fvMatrix<tensor>& m,
        const dictionary& dict
    ) const
    {
        return oversetFvMeshBase::solveOverset<tensor>(m, dict);
    }


    virtual bool update();

    virtual bool writeObject
    (
        IOstreamOption streamOpt,
        const bool writeOnProc
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
