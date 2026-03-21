# of-plugins

Standalone OpenFOAM plugin collection for rotating FSI cases.

This repository provides three independent plugins:

1. `solidBodyDisplacementLaplacianZone`
2. `dynamicOversetMotionSolverFvMesh`
3. `fsiOmega` (with `preciceOmega` Function1)

No submodules are used. The code is self-contained in this repository and built
as user libraries with `wmake`.

## Contents

- `solidBodyDisplacementLaplacianZone/`
  - Motion solver that combines:
    - rigid rotation of a selected zone (`cellZone` or `cellSet`), and
    - Laplacian mesh deformation for the rest of the mesh.
- `dynamicOversetMotionSolverFvMesh/`
  - Overset-enabled `dynamicFvMesh` wrapper for a single motion solver.
  - Lets you keep the same `motionSolver` and `...Coeffs` syntax used by
    `dynamicMotionSolverFvMesh`, but on overset meshes.
- `fsiOmega/`
  - `preciceOmega` Function1 that reads angular velocity from a
    `uniformDimensionedScalarField` (typically updated by preCICE adapter).

## Requirements

- OpenFOAM (tested with OpenFOAM.com line, v2406+ recommended)
- Standard OpenFOAM build environment loaded:
  - `WM_PROJECT_DIR`
  - `FOAM_USER_LIBBIN`
  - `wmake` in PATH

Optional for FSI coupling:
- preCICE adapter that writes angular velocity field (e.g. field name `omega`)

## Build

From repository root:

```bash
./Allwmake
```

This builds:

- `libsolidBodyDisplacementLaplacianZoneFvMotionSolver.so`
- `libdynamicOversetMotionSolverFvMesh.so`
- `libfsiOmega.so`

into `FOAM_USER_LIBBIN`.

You can also build each plugin separately:

```bash
cd solidBodyDisplacementLaplacianZone && ./Allwmake
cd ../dynamicOversetMotionSolverFvMesh && ./Allwmake
cd ../fsiOmega && ./Allwmake
```

## Usage: zonal rigid rotation + mesh deformation

Example `constant/dynamicMeshDict`:

```c++
dynamicFvMesh   dynamicMotionSolverFvMesh;

motionSolverLibs
(
    "libsolidBodyDisplacementLaplacianZoneFvMotionSolver.so"
    "libfvMotionSolvers.so"
);

motionSolver    solidBodyDisplacementLaplacianZone;

solidBodyDisplacementLaplacianZoneCoeffs
{
    // Use cellZone OR cellSet (not both)
    cellZone    rotorZone;
    // cellSet  rotorCells;

    solidBodyMotionFunction rotatingMotion;

    rotatingMotionCoeffs
    {
        origin  (0 0 0);
        axis    (0 1 0);
        omega   158;   // rad/s or any Function1<scalar>
    }

    diffusivity quadratic inverseDistance (rotorTip);
}
```

## Usage with preCICE angular velocity (FSI)

Load both libraries:

```c++
motionSolverLibs
(
    "libsolidBodyDisplacementLaplacianZoneFvMotionSolver.so"
    "libfvMotionSolvers.so"
    "libfsiOmega.so"
);
```

Then in `rotatingMotionCoeffs`:

```c++
omega   preciceOmega;

preciceOmegaCoeffs
{
    fieldName   omega;
}
```

Ensure your adapter configuration uses the same field name
(e.g. `nameOmegaField omega`).

## Usage on overset meshes

The important change for overset is the mesh class, not the motion-solver
syntax. In `constant/dynamicMeshDict` use:

```c++
dynamicFvMesh   dynamicOversetMotionSolverFvMesh;

motionSolverLibs
(
  "libsolidBodyDisplacementLaplacianZoneFvMotionSolver.so"
  "libfvMotionSolvers.so"
  "libfsiOmega.so"
);

motionSolver    solidBodyDisplacementLaplacianZone;

solidBodyDisplacementLaplacianZoneCoeffs
{
  cellZone    bladeZone;

  solidBodyMotionFunction rotatingMotion;

  rotatingMotionCoeffs
  {
    origin  (0 0 0);
    axis    (0 1 0);
    omega   158;

    // Or with FSI:
    // omega   preciceOmega;
    // preciceOmegaCoeffs
    // {
    //     fieldName omega;
    // }
  }

  diffusivity inverseDistance (propellerTip);
}
```

And in `system/controlDict` load the overset wrapper library in addition to the
standard overset support, for example:

```c++
libs (overset fvMotionSolvers dynamicOversetMotionSolverFvMesh);
```

This avoids the `solvers { ... }` multi-motion setup and keeps the same single
solver syntax as the non-overset plugin.

## Notes

- If neither `cellZone` nor `cellSet` is provided, rigid transform is applied
  to the full mesh.
- This repository intentionally contains only these plugins and no dependency
  on git submodules.

## License

GPL-3.0-or-later (same family as OpenFOAM user-extension code).
