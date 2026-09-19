# of-plugins

OpenFOAM plugins for **fluid–structure interaction (FSI)** of rotating
machinery: zonal rigid-body motion, Laplacian mesh deformation, overset
meshes, preCICE coupling, and actuator-line turbine loading.

![OpenFOAM](https://img.shields.io/badge/OpenFOAM-v2406%E2%80%93v2512-brightgreen)
![preCICE](https://img.shields.io/badge/preCICE-v3-blue)
![License](https://img.shields.io/badge/License-GPL--3.0--or--later-green)

## How the plugins fit together

1. **`precice-openfoam-adapter`** — preCICE coupling functionObject: exchanges
   forces and displacements with an external structural solver.
2. **`solidBodyDisplacementLaplacianZone`** — mesh-motion solver: applies the
   structural displacement to the blade surface and rotates the rotor zone
   rigidly, diffusing the deformation through the mesh with a Laplacian.
3. **`dynamicOversetZoneDisplacementFvMesh`** — runs the same motion solver on
   overset (chimera) meshes.
4. **`fsiOmega`** — `preciceOmega` Function1: rotor angular velocity read from
   preCICE for `rotatingMotion`.
5. **`turbinesFoam`** — actuator-line `fvOptions` adding aerodynamic loading to
   any compatible solver (`simpleFoam`, `pimpleFoam`, ...).

No submodules are used. All code is self-contained and built with `wmake`.

## Contents

- `solidBodyDisplacementLaplacianZone/`
  - Motion solver that combines:
    - rigid rotation of a selected zone (`cellZone` or `cellSet`), and
    - Laplacian mesh deformation for the rest of the mesh.
- `dynamicOversetZoneDisplacementFvMesh/`
  - Overset-enabled `dynamicFvMesh` wrapper for a single motion solver.
  - Lets you keep the same `motionSolver` and `...Coeffs` syntax used by
    `dynamicMotionSolverFvMesh`, but on overset meshes.
- `fsiOmega/`
  - `preciceOmega` Function1 that reads angular velocity from a
    `uniformDimensionedScalarField` (typically updated by preCICE adapter).
- `precice-openfoam-adapter/`
  - **Diverging fork** of [precice/openfoam-adapter](https://github.com/precice/openfoam-adapter)
  - **FSI-only version**: CHT (heat transfer) and FF (free-surface) modules removed.
  - Embedded here for simultaneous development with custom FSI modifications.
  - Includes **v2506 fix**: deprecated `Pstream::scatterList`/`gatherList` replaced with `OPstream`/`IPstream`.
  - See [precice-openfoam-adapter/README.md](precice-openfoam-adapter/README.md) for original documentation.
  - **Note**: This copy is independent — changes here do **not** automatically sync upstream.
- `turbinesFoam/`
  - **Vendored independent fork** of [turbinesFoam/turbinesFoam](https://github.com/turbinesFoam/turbinesFoam)
    (actuator line method for wind and marine turbines).
  - Original git history stripped; included as a base for research development.
  - Upstream `foamStyleCheck` tool removed (community vera++ checker, not an
    official OpenFOAM formatter) — the repository remains submodule-free.
  - **Fork extension:** optional blade **actuator surface model (ASM)**
    (`elementType actuatorSurfaceElement;` + `nChordwise` per line/blade)
    with chord-averaged inflow and mesh-based projection width; ALM remains
    the byte-identical default. HAWT tutorial at
    `tutorials/axialFlowTurbineASM/` with an ALM-vs-ASM comparison script.
  - See [turbinesFoam/README.md](turbinesFoam/README.md) for original documentation.

## FSI Motion Model

The `solidBodyDisplacementLaplacianZone` solver implements the correct FSI motion
model for rotating machinery:

```
U_total(x) = R(t) · (x + u_fsi(x)) - x
```

Where:
- `x` — original point coordinates
- `u_fsi` — structural deformation from FSI solver (applied via `cellDisplacement` BC)
- `R(t)` — rigid rotation matrix (from `solidBodyMotionFunction`, e.g. `rotatingMotion`)

**The order of operations is critical:**
1. **Deform first:** `x' = x + u_fsi` — apply structural deformation
2. **Rotate second:** `x_final = R(t) · x'` — apply rigid rotation to deformed body

This matches the physical behavior of a rotating blade: the blade deflects under
aerodynamic loads, then the whole (now-deformed) blade rotates.

The Laplacian diffusion (`div(gamma nabla U) = 0`) propagates `u_fsi` from the reference
patch (specified in `inverseDistance(patchName)`) through the rotation zone.
Cells at the zone boundary are pinned to zero displacement to prevent diffusion
outside the zone.

### AMI boundary protection: `boundaryDecay` diffusivity

When using AMI (Arbitrary Mesh Interface) meshes, non-linear diffusivity
functions can propagate deformation to the AMI boundary, breaking interface
connectivity. The `boundaryDecay` manipulator prevents this by smoothly
decaying the diffusivity to zero near selected patches.

```
patch AMI        decayDistance          interior
   |                 |
   |  D=0  --smooth--  D=D_base
   |                 |
```

Use it by wrapping your existing diffusivity chain:

```c++
diffusivity  boundaryDecay 0.05 (AMI.*) quadratic inverseDistance (blade);
```

| Parameter | Example | Description |
|-----------|---------|-------------|
| decayDistance | `0.05` | Distance (m) from AMI where D ramps 0 to 1 |
| patchNames | `(AMI.*)` | Regex matching AMI boundary patches |
| baseDiffusivity | `quadratic inverseDistance (blade)` | Any existing motionDiffusivity chain |

The decay factor uses a quintic smooth-step (`6x^5 - 15x^4 + 10x^3`) to ensure
C2 continuity (zero first and second derivatives) at both ends of the transition.

### AMI displacement clamping: `displacementDecay`

`boundaryDecay` prevents deformation from *reaching* the AMI, but the
cyclicAMI boundary condition can still interpolate non-zero displacement
back onto AMI points during `correctBoundaryConditions()`.  The optional
`displacementDecay` sub-dictionary clamps the solved displacement as a
post-processing step, guaranteeing zero motion at the interface:

```
patch AMI          decayDistance          interior
   |                   |
   |  U=0  --quintic--  U=U_solved
   |                   |
```

```c++
displacementDecay
{
    distance  0.05;      // decay length (m)
    patches   (AMI.*);   // regex matching AMI boundary patches
}
```

| Parameter | Example | Description |
|-----------|---------|-------------|
| distance | `0.05` | Distance (m) from the patches where displacement ramps 0 → 1 |
| patches | `(AMI.*)` | Regex matching the patches near which to clamp |

Behaviour:
- Points/cells **outside** the rotation zone: displacement set to exactly zero
- Points/cells **inside** the zone: quintic smooth-step (`6x^5 - 15x^4 + 10x^3`,
  C2 continuous) decay over `distance` measured from the patches via `patchWave`
- Applied after `correctBoundaryConditions()` in `curPoints()` and after the
  solve in `solve()`, so cyclicAMI interpolation cannot reintroduce non-zero
  displacement at the interface
- Skipped when no zone is selected (`moveAllCells` mode)

## Requirements

- OpenFOAM (tested with OpenFOAM.com line, v2406–v2512)
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

This builds all five plugins:

- `libsolidBodyDisplacementLaplacianZoneFvMotionSolver.so`
- `libdynamicOversetZoneDisplacementFvMesh.so`
- `libfsiOmega.so`
- `libpreciceAdapterFunctionObject.so` (adapter — warns but continues if
  preCICE's `pkg-config` file is missing)
- `libturbinesFoam.so`

into `FOAM_USER_LIBBIN`.

### Building the preCICE adapter (optional)

To build only the preCICE adapter (which has its own build system):

```bash
cd precice-openfoam-adapter && ./Allwmake
```

You can also build each component separately:

```bash
cd solidBodyDisplacementLaplacianZone && ./Allwmake
cd ../dynamicOversetZoneDisplacementFvMesh && ./Allwmake
cd ../fsiOmega && ./Allwmake
cd ../precice-openfoam-adapter && ./Allwmake   # requires preCICE dev files
cd ../turbinesFoam && ./Allwmake
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

    // For AMI meshes, wrap with boundaryDecay to protect the interface:
    // diffusivity boundaryDecay 0.05 (AMI.*) quadratic inverseDistance (rotorTip);

    // Optionally clamp displacement near the AMI after the solve:
    // displacementDecay
    // {
    //     distance 0.05;
    //     patches (AMI.*);
    // }
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

Same as above, but in a `dynamicOversetZoneDisplacementFvMesh` configuration.
In `controlDict` add:

```c++
libs (overset fvMotionSolvers dynamicOversetZoneDisplacementFvMesh);
```

Ensure your adapter configuration uses the same field name
(e.g. `nameOmegaField omega` in your `preciceDict`).

## Usage on overset meshes (different configuration)

Use this configuration for **overset (multi-mesh)** cases. The key differences are:

1. **dynamicFvMesh type** changes to `dynamicOversetZoneDisplacementFvMesh`
2. **controlDict** must load the overset wrapper library
3. **Otherwise**, the motion solver syntax remains identical

Example `constant/dynamicMeshDict`:

```c++
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

    // Optionally clamp displacement near the AMI after the solve:
    // displacementDecay
    // {
    //     distance 0.05;
    //     patches (AMI.*);
    // }
}
```

In `system/controlDict`, **must** load the overset wrapper library:

```c++
libs (overset fvMotionSolvers dynamicOversetZoneDisplacementFvMesh);
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
- `turbinesFoam` keeps its upstream script names: build with `./Allwmake`,
  clean with `./Allwclean` (the root `./Allclean` handles this automatically).

## License

GPL-3.0-or-later (same family as OpenFOAM user-extension code).
