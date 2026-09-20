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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::fv::nacelleSurfaceSampler::nNodes() const
{
    return positions_.size();
}


inline const Foam::List<Foam::point>&
Foam::fv::nacelleSurfaceSampler::positions() const
{
    return positionsBody_;
}


inline const Foam::List<Foam::vector>&
Foam::fv::nacelleSurfaceSampler::forces() const
{
    return forces_;
}


inline const Foam::List<Foam::vector>&
Foam::fv::nacelleSurfaceSampler::normals() const
{
    return normalsBody_;
}


inline const Foam::List<Foam::scalar>&
Foam::fv::nacelleSurfaceSampler::areas() const
{
    return areas_;
}


inline Foam::scalar Foam::fv::nacelleSurfaceSampler::rhoRef() const
{
    return rhoRef_;
}


inline Foam::scalar
Foam::fv::nacelleSurfaceSampler::referenceVelocity() const
{
    return referenceVelocity_;
}


inline const Foam::List<Foam::point>&
Foam::fv::nacelleSurfaceSampler::positionsGlobal() const
{
    return positions_;
}


inline const Foam::List<Foam::vector>&
Foam::fv::nacelleSurfaceSampler::normalsGlobal() const
{
    return normals_;
}


inline const Foam::tensor&
Foam::fv::nacelleSurfaceSampler::bodyToGlobal() const
{
    return bodyToGlobal_;
}


// ************************************************************************* //
