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

inline const Foam::List<Foam::point>&
Foam::fv::bladeSurfaceSampler::canonicalNodes() const
{
    return canonicalNodes_;
}


inline const Foam::List<Foam::scalar>&
Foam::fv::bladeSurfaceSampler::station() const
{
    return station_;
}


inline const Foam::List<Foam::scalar>&
Foam::fv::bladeSurfaceSampler::chordFraction() const
{
    return chordFraction_;
}


inline const Foam::List<Foam::label>&
Foam::fv::bladeSurfaceSampler::patch() const
{
    return patch_;
}


inline const Foam::List<Foam::scalar>&
Foam::fv::bladeSurfaceSampler::patchAreaShare() const
{
    return patchAreaShare_;
}


inline const Foam::List<Foam::scalar>&
Foam::fv::bladeSurfaceSampler::nodeH() const
{
    return nodeH_;
}


inline const Foam::List<Foam::scalar>&
Foam::fv::bladeSurfaceSampler::elementStation() const
{
    return elementStation_;
}


inline const Foam::List<Foam::scalar>&
Foam::fv::bladeSurfaceSampler::patchArea() const
{
    return patchArea_;
}


inline const Foam::List<Foam::vector>&
Foam::fv::bladeSurfaceSampler::forces() const
{
    return forces_;
}


inline const Foam::word& Foam::fv::bladeSurfaceSampler::kernelMode() const
{
    return kernelMode_;
}


inline Foam::scalar Foam::fv::bladeSurfaceSampler::meshFactor() const
{
    return meshFactor_;
}


// ************************************************************************* //
