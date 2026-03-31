# AGENTS.md — of-plugins

OpenFOAM plugin collection for rotating FSI cases and preCICE coupling. C++ codebase using OpenFOAM's `wmake` build system.

## Build / Clean Commands

### Build all plugins
```bash
./Allwmake                    # builds all 4 plugins
```

### Build a single plugin
```bash
cd solidBodyDisplacementLaplacianZone && ./Allwmake
cd dynamicOversetZoneDisplacementFvMesh && ./Allwmake
cd fsiOmega && ./Allwmake
cd precice-openfoam-adapter && ./Allwmake   # requires preCICE installed
```

### Clean
```bash
./Allclean                    # cleans all plugins
wclean                        # cleans current plugin directory
```

### Clean individual plugin
```bash
cd <plugin-dir> && ./Allclean
```

### Parallel builds
Set `WM_NCOMPPROCS` before building:
```bash
export WM_NCOMPPROCS=4
./Allwmake
```

### Build output
Libraries are installed to `$FOAM_USER_LIBBIN`.

## Testing

**No automated test suite exists.** This repository contains no unit tests or integration tests. Validation is done manually by running OpenFOAM cases that use the plugins.

When adding new code:
- Verify compilation with `./Allwmake`
- Test with a real OpenFOAM case if possible
- For MPI-related changes, test in parallel (`mpirun -np N`)

## Code Style

Follow **OpenFOAM coding conventions**. The `.clang-format` in `precice-openfoam-adapter/` is the authoritative style reference.

### Formatting
- **Indentation:** 4 spaces, no tabs
- **Braces:** Allman style (opening brace on new line)
- **Line length:** No hard limit (ColumnLimit: 0)
- **Max empty lines:** 2
- **Include sorting:** Disabled — preserve OpenFOAM include order

### File header
Every source file must start with the OpenFOAM license block:
```
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) <YEAR> <AUTHOR>
-------------------------------------------------------------------------------
License
    GPL-3.0-or-later
\*---------------------------------------------------------------------------*/
```

### Includes
- `.C` files: include the matching `.H` first, then OpenFOAM headers alphabetically
- Use OpenFOAM include paths: `-I$(LIB_SRC)/finiteVolume/lnInclude`, etc.
- Do NOT sort includes automatically (breaks OpenFOAM lnInclude ordering)

### Naming conventions
- **Classes:** `camelCase` starting with lowercase (e.g., `solidBodyDisplacementLaplacianZoneFvMotionSolver`)
- **Member variables:** `camelCase` with trailing underscore for private data (e.g., `cellDisplacement_`, `diffusivityPtr_`)
- **Functions/methods:** `camelCase` (e.g., `curPoints()`, `updateMesh()`)
- **Local variables:** `camelCase` (e.g., `cellIDs`, `zoneID`)
- **Constants:** `UPPER_CASE` or OpenFOAM constants via `Foam::constant::mathematical::pi`

### Types
- Use OpenFOAM types: `label` (not `int`), `scalar` (not `double`), `vector`, `tensor`
- Use `word` for string identifiers
- Use `tmp<Field>` for temporary field returns
- Use `autoPtr` / `refPtr` for owned pointers where appropriate

### Runtime type registration
```cpp
namespace Foam
{
    defineTypeNameAndDebug(ClassName, 0);
    addToRunTimeSelectionTable(baseClass, ClassName, dictionary);
}
```

### Error handling
- Use `FatalErrorInFunction` / `FatalIOErrorInFunction` for fatal errors
- Use `WarningInFunction` for non-fatal warnings
- Use `Info<< ... << endl;` for informational messages
- Use `returnReduce()` for parallel-safe reductions
- Check return values from `findZoneID()`, `findObject<>()`, etc.

### MPI / parallel considerations
- Use `OPstream`/`IPstream` (NOT deprecated `Pstream::scatterList`/`gatherList`)
- Call `syncTools::syncPointList()` / `syncTools::syncFaceList()` for parallel consistency
- Use `returnReduce()` with `sumOp<label>()` etc. for global reductions
- Be aware of overset mesh MPI communication patterns — they can corrupt heap if misused

### Overset mesh specifics
- Use `updateCoeffs()` (not `correctBoundaryConditions()`) on overset patches
- Defer diffusivity creation to avoid constructor-order crashes
- Zero `cellDisplacement_` after mesh motion to prevent checkpoint corruption
- Pin hole cells in matrix with `VGREAT` to prevent Laplacian singularity

### Mathematical constants
- Use `Foam::constant::mathematical::pi` (NOT `M_PI`)
- Include `mathematicalConstants.H` when needed

## Shell scripts
- All `Allwmake`/`Allclean` scripts use `#!/bin/sh` (POSIX sh, not bash)
- Use `cd "${0%/*}" || exit` pattern for directory change
- Use `set -e -u` in clean scripts

## Known issues
- Root `Allwmake`/`Allclean` reference `dynamicOversetMotionSolverFvMesh` but the directory is named `dynamicOversetZoneDisplacementFvMesh`. Fix these when touching the scripts.

## Git
- Do not commit build artifacts (`*.o`, `*.so`, `platforms/`, `lnInclude/`, `*.log`)
- License: GPL-3.0-or-later
