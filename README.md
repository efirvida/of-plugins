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

## Usage: zonal rigid rotation + mesh deformation (non-overset)

Use this configuration for **standard (single-mesh, no overset)** cases.

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

In `system/controlDict` add standard motion solver libraries:

```c++
libs (fvMotionSolvers);
```

## Usage with preCICE angular velocity (FSI coupling)

The `preciceOmega` Function1 works with **both** non-overset and overset configurations.

### With non-overset meshes:

Load both libraries in `motionSolverLibs`:

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

### With overset meshes:

Same as above, but in a `dynamicOversetMotionSolverFvMesh` configuration.
In `controlDict` add:

```c++
libs (overset fvMotionSolvers dynamicOversetMotionSolverFvMesh);
```

Ensure your adapter configuration uses the same field name
(e.g. `nameOmegaField omega` in your `preciceDict`).

## Usage on overset meshes (different configuration)

Use this configuration for **overset (multi-mesh)** cases. The key differences are:

1. **dynamicFvMesh type** changes to `dynamicOversetMotionSolverFvMesh`
2. **controlDict** must load the overset wrapper library
3. **Otherwise**, the motion solver syntax remains identical

Example `constant/dynamicMeshDict`:

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

In `system/controlDict`, **must** load the overset wrapper library:

```c++
libs (overset fvMotionSolvers dynamicOversetMotionSolverFvMesh);
```

**Key point:** This avoids the `solvers { ... }` multi-motion accumulation approach that was problematic before. The overset wrapper allows you to keep **a single motion solver** with simple syntax, even on overset meshes.

## Summary of configuration changes

| Aspect | Non-overset | Overset |
|--------|-------------|---------|
| `dynamicFvMesh` | `dynamicMotionSolverFvMesh` | `dynamicOversetMotionSolverFvMesh` |
| `motionSolver` | `solidBodyDisplacementLaplacianZone` | `solidBodyDisplacementLaplacianZone` (same) |
| `...Coeffs` | `solidBodyDisplacementLaplacianZoneCoeffs` | `solidBodyDisplacementLaplacianZoneCoeffs` (same) |
| controlDict libs | `(fvMotionSolvers)` | `(overset fvMotionSolvers dynamicOversetMotionSolverFvMesh)` |
| preCICE support | Yes, add `libfsiOmega` | Yes, same as non-overset |

## Notes

- If neither `cellZone` nor `cellSet` is provided, rigid transform is applied
  to the full mesh (both non-overset and overset).
- The overset wrapper eliminates the need for `solvers { ... }` dictionaries
  and the problematic multi-solver accumulation approach.
- This repository intentionally contains only these plugins and has no dependency
  on git submodules.

## License

GPL-3.0-or-later (same family as OpenFOAM user-extension code).
