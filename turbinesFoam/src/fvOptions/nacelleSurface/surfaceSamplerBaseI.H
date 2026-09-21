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

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

inline Foam::label Foam::fv::surfaceSamplerBase::nNodes() const
{
    return positions_.size();
}


inline const Foam::List<Foam::point>&
Foam::fv::surfaceSamplerBase::positions() const
{
    return positionsBody_;
}


inline const Foam::List<Foam::vector>&
Foam::fv::surfaceSamplerBase::normals() const
{
    return normalsBody_;
}


inline const Foam::List<Foam::scalar>&
Foam::fv::surfaceSamplerBase::areas() const
{
    return areas_;
}


inline const Foam::List<Foam::point>&
Foam::fv::surfaceSamplerBase::positionsGlobal() const
{
    return positions_;
}


inline const Foam::List<Foam::vector>&
Foam::fv::surfaceSamplerBase::normalsGlobal() const
{
    return normals_;
}


inline const Foam::tensor&
Foam::fv::surfaceSamplerBase::bodyToGlobal() const
{
    return bodyToGlobal_;
}


inline Foam::scalar Foam::fv::surfaceSamplerBase::rhoRef() const
{
    return rhoRef_;
}


// ************************************************************************* //
