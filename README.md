# of-plugins

Standalone OpenFOAM plugin collection for rotating FSI cases and preCICE coupling.

This repository provides four independent plugins:

1. `solidBodyDisplacementLaplacianZone` — motion solver for zonal rigid body + mesh deformation
2. `dynamicOversetMotionSolverFvMesh` — overset wrapper for the above
3. `fsiOmega` — `preciceOmega` Function1 for preCICE angular velocity coupling
4. `precice-openfoam-adapter` — preCICE OpenFOAM adapter (diverging fork)

No submodules are used. All code is self-contained and built with `wmake`.

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
- `precice-openfoam-adapter/`
  - **Diverging fork** of [precice/openfoam-adapter](https://github.com/precice/openfoam-adapter)
  - Embedded here for simultaneous development and modifications.
  - See [precice-openfoam-adapter/README.md](precice-openfoam-adapter/README.md) for original documentation.
  - **Note**: This copy is independent — changes here do **not** automatically sync upstream.

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

This builds the three motion-related libraries:

- `libsolidBodyDisplacementLaplacianZoneFvMotionSolver.so`
- `libdynamicOversetMotionSolverFvMesh.so`
- `libfsiOmega.so`

into `FOAM_USER_LIBBIN`.

### Building the preCICE adapter (optional)

To build only the preCICE adapter (which has its own build system):

```bash
cd precice-openfoam-adapter && ./Allwmake
```

You can also build each component separately:

```bash
cd solidBodyDisplacementLaplacianZone && ./Allwmake
cd ../dynamicOversetMotionSolverFvMesh && ./Allwmake
cd ../fsiOmega && ./Allwmake
cd ../precice-openfoam-adapter && ./Allwmake   # requires preCICE dev files
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

## About the preCICE adapter divergence

The `precice-openfoam-adapter/` subdirectory is a **complete independent copy** of the 
[precice/openfoam-adapter](https://github.com/precice/openfoam-adapter) repository
at a point-in-time snapshot. It is **not** a git submodule.

**Why?** To allow simultaneous, independent development of custom modifications (e.g., enhanced
FSI control, overset integration, tighter coupling with `preciceOmega`) without depending on
upstream releases.

**Important:**
- Changes made here **will diverge** from the upstream precice/openfoam-adapter main branch.
- To sync back with upstream or integrate upstream changes, you must manually handle merges.
- For bug fixes or features that should go upstream, consider opening PRs on the
  [precice/openfoam-adapter](https://github.com/precice/openfoam-adapter) repo directly.

## Notes

- If neither `cellZone` nor `cellSet` is provided, rigid transform is applied
  to the full mesh (both non-overset and overset).
- The overset wrapper eliminates the need for `solvers { ... }` dictionaries
  and the problematic multi-solver accumulation approach.
- This repository intentionally contains only these plugins and has no dependency
  on git submodules.

## License

GPL-3.0-or-later (same family as OpenFOAM user-extension code).
