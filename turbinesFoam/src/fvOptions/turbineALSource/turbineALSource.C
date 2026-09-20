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

#include "turbineALSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"

#include <cmath>

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(turbineALSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        turbineALSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::turbineALSource::rotateVector
(
    vector& vectorToRotate,
    vector rotationPoint,
    vector axis,
    scalar radians
)
{
    // Declare and define the rotation matrix (from SOWFA)
    tensor RM;
    scalar angle = radians;
    RM.xx() = Foam::sqr(axis.x())
            + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle);
    RM.xy() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle);
    RM.xz() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle);
    RM.yy() = Foam::sqr(axis.y())
            + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z())
            + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    // Rotation matrices make a rotation about the origin, so need to subtract
    // rotation point off the point to be rotated.
    vectorToRotate -= rotationPoint;

    // Perform the rotation.
    vectorToRotate = RM & vectorToRotate;

    // Return the rotated point to its new location relative to the rotation
    // point
    vectorToRotate += rotationPoint;
}


void Foam::fv::turbineALSource::createCoordinateSystem()
{
    // Should be unique for each type of turbine
}


void Foam::fv::turbineALSource::createBlades()
{
    // Should be unique for each type of turbine
}


void Foam::fv::turbineALSource::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = time_.path()/"../postProcessing/turbines"
            / time_.timeName();
    }
    else
    {
        dir = time_.path()/"postProcessing/turbines"
            / time_.timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_ = new OFstream(dir/name_ + ".csv");

    *outputFile_<< "time,angle_deg,tsr,cp,cd,ct";

    forAll(blades_, i)
    {
        *outputFile_<< ",cd_" << bladeNames_[i];
        *outputFile_<< ",ct_" << bladeNames_[i];
    }

    *outputFile_<< endl;
}


// W3.2 registry conventions, resolved against the seams already in this tree:
// the registry host is Time, matching the adapter's global-data coupling
// (modules/generic/Generic.C uses mesh_.time().foundObject/lookupObject and
// sortedNames<uniformDimensionedScalarField>(), ReadWrite.C uses
// runTime_.foundObject/lookupObjectRef) and fsiOmega/preciceOmega.C
// (runTime.foundObject/lookupObject). Those seams only ever create the field
// in memory, so they give no precedent for the optional disk read-back: the
// OpenFOAM idiom for "read it if the start time directory has it, otherwise
// default" is IOobject::READ_IF_PRESENT at the current time instance, which
// under startFrom latestTime is the latest written time directory.
void Foam::fv::turbineALSource::createAngleDegField()
{
    const word fieldName("angleDeg." + name_);

    IOobject io
    (
        fieldName,
        time_.timeName(),
        time_,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    // Restart/rollback: a field written at the start time is re-read from
    // that time directory (READ_IF_PRESENT)
    const bool haveAngle = io.typeHeaderOk<uniformDimensionedScalarField>(true);

    angleDegField_.reset
    (
        new uniformDimensionedScalarField
        (
            io,
            dimensionedScalar(fieldName, dimless, scalar(0))
        )
    );

    if (haveAngle)
    {
        angle0_ = angleDegField_->value();
        t0_ = time_.value();

        Info<< "Resuming azimuth of " << name_ << " from " << fieldName
            << ": " << angle0_ << " deg at t = " << t0_ << endl;
    }
    else
    {
        angle0_ = 0.0;
        t0_ = time_.startTime().value();
    }

    // The last applied azimuth initialises to the interval start
    angleDeg_ = angle0_;
}


void Foam::fv::turbineALSource::createOmegaOverrideField()
{
    if (time_.foundObject<uniformDimensionedScalarField>(omegaOverrideField_))
    {
        return;
    }

    // Register in the Time registry so the preCICE adapter's globalData
    // ReadWrite can bind to it by name (the field is the seam; the adapter
    // creates/owns its own copy when this source does not)
    omegaOverrideFieldPtr_.reset
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                omegaOverrideField_,
                time_.constant(),
                time_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dimensionedScalar
            (
                omegaOverrideField_,
                dimensionSet(0, 0, -1, 0, 0, 0, 0),
                scalar(0)
            )
        )
    );

    Info<< "Created omega override field " << omegaOverrideField_ << endl;
}


void Foam::fv::turbineALSource::updateTSROmega()
{
    if (hasOmegaOverride_)
    {
        // Angular velocity [rad/s] supplied through the registry override
        // field (the preCICE adapter's globalData seam, mirroring fsiOmega)
        omega_ = time_.lookupObject<uniformDimensionedScalarField>
        (
            omegaOverrideField_
        ).value();
        tipSpeedRatio_ = omega_*rotorRadius_/mag(freeStreamVelocity_);
    }
    else
    {
        // Update tip speed ratio and omega
        scalar theta = degToRad(angleDeg_);
        tipSpeedRatio_ = meanTSR_ + tsrAmplitude_*cos(nBlades_*(theta - tsrPhase_));
        omega_ = tipSpeedRatio_*mag(freeStreamVelocity_)/rotorRadius_;
    }
}


Foam::scalar Foam::fv::turbineALSource::azimuth(const scalar t) const
{
    const scalar pi = constant::mathematical::pi;
    const scalar omega0 = meanTSR_*mag(freeStreamVelocity_)/rotorRadius_;
    const scalar theta0 = degToRad(angle0_);

    // Constant TSR: theta(t) = theta0 + omega0*(t - t0)
    scalar theta = theta0 + omega0*(t - t0_);

    if (tsrAmplitude_ != 0.0)
    {
        // Oscillating TSR: dtheta/dt = omega0*(1 + m*cos(nB*(theta - phi)))
        // is separable. With u = nB*(theta - phi):
        //
        //     G(u/2) = G(u0/2) + (nB*sqrt(1 - m^2)/2)*omega0*(t - t0)
        //
        // where G(x) = atan(k*tan(x)) + pi*round(x/pi) is the continuous
        // branch of the antiderivative and k = sqrt((1 - m)/(1 + m))
        if (mag(tsrAmplitude_) >= mag(meanTSR_) || mag(meanTSR_) <= VSMALL)
        {
            FatalErrorInFunction
                << "The closed-form tsrAmplitude oscillation requires "
                << "|tsrAmplitude| < |tipSpeedRatio|" << nl
                << "    tipSpeedRatio: " << meanTSR_ << nl
                << "    tsrAmplitude:  " << tsrAmplitude_
                << exit(FatalError);
        }

        const scalar m = tsrAmplitude_/meanTSR_;
        const scalar k = Foam::sqrt((1.0 - m)/(1.0 + m));
        const scalar nB = scalar(nBlades_);

        const scalar u0 = nB*(theta0 - tsrPhase_);
        const scalar x0 = 0.5*u0;
        const scalar y0 = Foam::atan(k*Foam::tan(x0)) + pi*std::round(x0/pi);

        const scalar y = y0 + 0.5*nB*Foam::sqrt(1.0 - sqr(m))*omega0*(t - t0_);

        // Invert G; round(y/pi) unwraps the atan branch so that theta(t)
        // stays monotonic and continuous
        const scalar n = std::round(y/pi);
        const scalar x = n*pi + Foam::atan(Foam::tan(y)/k);

        theta = tsrPhase_ + 2.0*x/nB;
    }

    // The azimuth is returned unwrapped (continuous and monotonic). It is the
    // same angle modulo 360 degrees, and keeping the integral single-valued is
    // what makes the CSV angle_deg column reproduce the accumulator beyond a
    // full revolution (spec: CSV output identical to the accumulator).
    return radToDeg(theta);
}


void Foam::fv::turbineALSource::rotate()
{
    const scalar t = time_.value();
    const scalar anglePrev = angleDeg_;

    if (hasOmegaOverride_)
    {
        // The override is an external stepwise signal: integrate it
        // incrementally, as its history cannot be recovered from
        // (t0_, angle0_, t). The interval is the time elapsed since the last
        // applied rotation, matching the accumulator this replaces (and
        // making a repeated call at the same time a no-op)
        angleDeg_ += radToDeg(omega_*(t - lastRotationTime_));
    }
    else
    {
        // Idempotent, time-derived azimuth: pure function of (t0_, angle0_, t)
        angleDeg_ = azimuth(t);
    }

    // Rotation to apply is the change since the last applied azimuth
    const scalar deltaRad = degToRad(angleDeg_ - anglePrev);
    rotate(deltaRad);

    lastRotationTime_ = t;
    updateTSROmega();

    // Persist every step (not only at write intervals) so a restart or a
    // preCICE rollback to an arbitrary time can restore the azimuth; these
    // per-step angleDeg.<name> files are the checkpoint record
    if (angleDegField_.valid())
    {
        angleDegField_->value() = angleDeg_;
        angleDegField_->write();
    }
}


void Foam::fv::turbineALSource::rotate(scalar radians)
{
    // Should be defined for each turbine type
}


void Foam::fv::turbineALSource::printPerf()
{
    Info<< "Azimuthal angle (degrees) of " << name_ << ": " << angleDeg_
        << endl;
    Info<< "Tip speed ratio of " << name_ << ": " << tipSpeedRatio_ << endl;
    Info<< "Power coefficient from " << name_ << ": " << powerCoefficient_
        << endl;
    Info<< "Rotor drag coefficient from " << name_ << ": " << dragCoefficient_
        << endl << endl;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::turbineALSource::turbineALSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    time_(mesh.time()),
    lastRotationTime_(time_.value()),
    rhoRef_(1.0),
    omega_(0.0),
    angleDeg_(0.0),
    angle0_(0.0),
    t0_(time_.startTime().value()),
    nBlades_(0),
    freeStreamVelocity_(vector::zero),
    forceField_
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "force",
            dimForce/dimVolume/dimDensity,
            vector::zero
        )
    ),
    frontalArea_(0.0),
    powerCoefficient_(0.0),
    dragCoefficient_(0.0),
    torqueCoefficient_(0.0),
    hasOmegaOverride_(false),
    omegaOverrideField_(word::null),
    omegaOverrideFieldPtr_(),
    angleDegField_()
{
    // Register the persisted azimuth field and seed the restart state
    createAngleDegField();

    forceField_.write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::turbineALSource::~turbineALSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::turbineALSource::angleDeg() const
{
    return angleDeg_;
}


void Foam::fv::turbineALSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Should be unique for each turbine type
}


void Foam::fv::turbineALSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Should be unique for each turbine type
}


void Foam::fv::turbineALSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // Should be unique for each turbine type
}


void Foam::fv::turbineALSource::printCoeffs() const
{
    Info<< "Number of blades: " << nBlades_ << endl;
}


void Foam::fv::turbineALSource::writePerf()
{
    *outputFile_<< time_.value() << "," << angleDeg_ << ","
                << tipSpeedRatio_ << "," << powerCoefficient_ << ","
                << dragCoefficient_ << "," << torqueCoefficient_;

    // Write power, drag, and torque coefficients for each blade
    forAll(blades_, i)
    {
        // Write drag (thrust) coefficient contribution from blade
        scalar bladeCd = blades_[i].force() & freeStreamDirection_
            / (0.5*frontalArea_*magSqr(freeStreamVelocity_));
        *outputFile_<< "," << bladeCd;
        // Write torque coefficient contribution from blade
        scalar bladeTorque = bladeMoments_[i] & axis_;
        scalar bladeCt = bladeTorque
            / (0.5*frontalArea_*rotorRadius_* magSqr(freeStreamVelocity_));
        *outputFile_<< "," << bladeCt;
    }

    *outputFile_<< endl;
}


void Foam::fv::turbineALSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::turbineALSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Read coordinate system/geometry invariant properties
        coeffs_.lookup("origin") >> origin_;
        coeffs_.lookup("axis") >> axis_;
        axis_ /= mag(axis_);
        coeffs_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
        coeffs_.lookup("tipSpeedRatio") >> meanTSR_;
        coeffs_.lookup("rotorRadius") >> rotorRadius_;
        tsrAmplitude_ = coeffs_.lookupOrDefault("tsrAmplitude", 0.0);
        tsrPhase_ = coeffs_.lookupOrDefault("tsrPhase", 0.0);

        // Optional registry omega override field (FSI seam)
        omegaOverrideField_ = coeffs_.getOrDefault<word>
        (
            "omegaOverrideField",
            word::null
        );
        hasOmegaOverride_ = !omegaOverrideField_.empty();
        if (hasOmegaOverride_)
        {
            createOmegaOverrideField();
        }

        // Get blade information
        bladesDict_ = coeffs_.subDict("blades");
        nBlades_ = bladesDict_.keys().size();
        bladeNames_ = bladesDict_.toc();
        bladeMoments_.setSize(nBlades_);

        // Set tip speed ratio and omega
        updateTSROmega();

        // Get dynamic stall subdict
        dynamicStallDict_ = coeffs_.subOrEmptyDict("dynamicStall");

        // Get profiles information
        profileData_ = coeffs_.subDict("profileData");

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
