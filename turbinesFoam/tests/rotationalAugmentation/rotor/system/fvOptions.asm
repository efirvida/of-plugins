/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2506                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvOptions;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Same augmentation-on rotor as `fvOptions.on`, but the elements are the
// no-mesh actuator surface element: the correction must be inherited through
// the shared force chain with no surface-specific code.

turbine
{
    type            axialFlowTurbineALSource;
    active          on;

    axialFlowTurbineALSourceCoeffs
    {
        fieldNames          (U);
        selectionMode       cellSet;
        cellSet             turbine;
        origin              (0 0 0);
        axis                (1 0 0);
        verticalDirection   (0 0 1);
        freeStreamVelocity  (10 0 0);
        tipSpeedRatio       3.0;
        rotorRadius         0.45;

        dynamicStall
        {
            active          off;
            dynamicStallModel LeishmanBeddoes;
        }

        rotationalAugmentation
        {
            active          on;
            model           DuSelig;
            a               1;
            b               1;
            d               1;
        }

        endEffects
        {
            active          off;
            endEffectsModel Glauert;
        }

        blades
        {
            blade1
            {
                writePerf           true;
                writeElementPerf    true;
                elementType         actuatorSurfaceElement;
                nChordwise          5;
                nElements           2;
                elementProfiles     (S809 S809);
                elementData
                (
                    (0.0  0.05  0.0  0.04  0.25  0.0)
                    (0.0  0.45  0.0  0.04  0.25  0.0)
                );
            }
        }

        profileData
        {
            S809
            {
                data (#include "S809.dat");
            }
        }
    }
}

// ************************************************************************* //
