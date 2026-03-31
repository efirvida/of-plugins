/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2026 OpenFOAM Foundation
    Copyright (C) 2024-2026 Contributors to the fsiOmega library
-------------------------------------------------------------------------------
License
    This file is part of the fsiOmega library for OpenFOAM.

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

#include "preciceOmega.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{
    defineTypeNameAndDebug(preciceOmega, 0);

    Function1<scalar>::adddictionaryConstructorToTable
        <FieldFunction1<preciceOmega>>
        addpreciceOmegascalarConstructorToTable_;
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::Function1Types::preciceOmega::read(const dictionary& dict)
{
    fieldName_ = dict.getOrDefault<word>("fieldName", "omega");
    
    // Reset field pointers since we might have a new field name
    omegaField_ = nullptr;
    omegaFieldOwning_.reset();

    if (debug)
    {
        Info<< "preciceOmega:" << nl
            << "    fieldName = " << fieldName_ << nl
            << "    Reading omega from uniformDimensionedScalarField" << nl
            << "    (updated by preCICE adapter)" << endl;
    }
}


void Foam::Function1Types::preciceOmega::initField() const
{
    if (omegaField_)
    {
        return;  // Already initialized
    }

    if (!obrPtr_)
    {
        FatalErrorInFunction
            << "preciceOmega requires access to the object registry" << nl
            << "    to read the omega field '" << fieldName_ << "'" << nl
            << "    Make sure the Function1 is constructed with a valid" << nl
            << "    objectRegistry pointer (e.g., from mesh or Time)."
            << abort(FatalError);
    }

    // The preCICE adapter registers the omega field in Time's objectRegistry.
    // We always search in Time to ensure we find the same field.
    const Time& runTime = obrPtr_->time();
    
    // Look for the field in Time registry (where the adapter creates it)
    if (runTime.foundObject<uniformDimensionedScalarField>(fieldName_))
    {
        omegaField_ = 
            &const_cast<uniformDimensionedScalarField&>(
                runTime.lookupObject<uniformDimensionedScalarField>(fieldName_));
        
        if (debug)
        {
            Info<< "preciceOmega: Found omega field '" << fieldName_ 
                << "' in Time registry" << endl;
        }
        return;
    }
    
    // Field not found - create it in Time registry
    // The adapter will find and use this same field when it initializes
    omegaFieldOwning_.reset(new uniformDimensionedScalarField(
        IOobject(
            fieldName_,
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(
            fieldName_,
            dimensionSet(0, 0, -1, 0, 0, 0, 0),
            0.0
        )
    ));

    omegaField_ = omegaFieldOwning_.get();
    
    if (debug)
    {
        Info<< "preciceOmega: Created omega field '" << fieldName_ 
            << "' in Time registry (adapter will use same field)" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1Types::preciceOmega::preciceOmega
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<scalar>(entryName, dict, obrPtr),
    obrPtr_(obrPtr),
    fieldName_("omega"),
    omegaField_(nullptr)
{
    read(dict);
}


Foam::Function1Types::preciceOmega::preciceOmega
(
    const preciceOmega& rhs
)
:
    Function1<scalar>(rhs),
    obrPtr_(rhs.obrPtr_),
    fieldName_(rhs.fieldName_),
    omegaField_(nullptr)  // Don't copy the pointer, will be looked up again
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1Types::preciceOmega::value(const scalar t) const
{
    initField();
    
    const scalar omega = omegaField_->value();
    const scalar rpm = omega * 60.0 / (2.0 * constant::mathematical::pi);

    Info<< "preciceOmega::value(t=" << t << ")"
        << " omega=" << omega << " [rad/s]"
        << " (" << rpm << " [rpm])" << endl;

    return omega;
}


Foam::scalar Foam::Function1Types::preciceOmega::integrate
(
    const scalar t1,
    const scalar t2
) const
{
    initField();
    
    const scalar omega = omegaField_->value();
    const scalar rpm = omega * 60.0 / (2.0 * constant::mathematical::pi);
    const scalar integral = omega * (t2 - t1);
    
    Info<< "preciceOmega::integrate(t1=" << t1 << ", t2=" << t2 << ")" << nl
        << "    omega=" << omega << " [rad/s] (" << rpm << " [rpm])" << nl
        << "    integral=" << integral << " [rad]" << endl;

    return integral;
}


void Foam::Function1Types::preciceOmega::writeData(Ostream& os) const
{
    Function1<scalar>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));
    os.writeEntry("fieldName", fieldName_);
    os.endBlock();
}


// ************************************************************************* //
