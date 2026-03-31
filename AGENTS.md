# AGENTS.md — of-plugins

OpenFOAM plugin collection for rotating FSI cases and preCICE coupling.
C++ codebase using OpenFOAM's `wmake` build system. Each plugin is an
independent shared library installed to `$FOAM_USER_LIBBIN`.

## Repository layout

```
solidBodyDisplacementLaplacianZone/   # Laplacian motion solver with rigid zone support
dynamicOversetZoneDisplacementFvMesh/ # Overset mesh wrapper for zone-based motion
fsiOmega/                             # Angular-velocity Function1 driven by preCICE
precice-openfoam-adapter/             # Full OpenFOAM–preCICE coupling adapter (FSI-only fork)
```

## Build / Clean

```sh
./Allwmake                             # build all plugins (serial)
WM_NCOMPPROCS=4 ./Allwmake             # build all plugins (parallel)
./Allclean                             # clean all plugins

# Build or clean a single plugin:
cd solidBodyDisplacementLaplacianZone && ./Allwmake
cd dynamicOversetZoneDisplacementFvMesh && ./Allwmake
cd fsiOmega && ./Allwmake
cd precice-openfoam-adapter && ./Allwmake   # requires preCICE + pkg-config
```

The adapter's `Allwmake` uses `pkg-config libprecice` automatically and emits
warnings (not errors) when the `.pc` file is missing. Environment variables:
- `PRECICE_OPENFOAM_CFLAGS` — extra preprocessor flags
- `PRECICE_OPENFOAM_TARGET_DIR` — install destination (default `$FOAM_USER_LIBBIN`)

## Testing

**No automated test suite.** Validation is manual — run a real OpenFOAM case
that loads the compiled library. When adding new code:
- `./Allwmake` must exit 0 with zero `error:` lines in `wmake.log`
- The adapter's `Allwmake` runs `ldd -r` and checks for undefined symbols
- For MPI-related changes, test in parallel with `mpirun -np N`

## CI/CD (precice-openfoam-adapter only)

Workflows live in `precice-openfoam-adapter/.github/workflows/`:
- `build.yml` — builds against OpenFOAM v2512 on push/PR
- `pre-commit.yml` — runs pre-commit hooks on push to main/develop and all PRs
- `check-format.yml` — runs clang-format v11 on all `.C`/`.H` files
- `check-shell.yml` — runs shellcheck on all shell scripts

## Code formatting

The authoritative style is `precice-openfoam-adapter/.clang-format` (clang-format 11).
```sh
cd precice-openfoam-adapter
sh tools/format-code.sh               # format all .C/.H files
```

Pre-commit hooks (clang-format, shellcheck, markdownlint) are configured in
`precice-openfoam-adapter/.pre-commit-config.yaml`:
```sh
pip install pre-commit && pre-commit install
```

## Code style

Follow **OpenFOAM coding conventions** throughout all four plugins.

### File header
Every `.C` and `.H` file must begin with the standard OpenFOAM license block
(GPL-3.0-or-later) including copyright year and author.

### Formatting
- **Indentation:** 4 spaces, no tabs
- **Braces:** Allman style — opening brace on its own line, always
- **Line length:** no hard limit (`ColumnLimit: 0`)
- **Blank lines:** at most 2 consecutive
- **Pointer/reference alignment:** left (`scalar* ptr`, not `scalar *ptr`)

### Includes
- In `.C` files: include the matching `.H` first, then other headers
- Use OpenFOAM lnInclude paths: `-I$(LIB_SRC)/finiteVolume/lnInclude`, etc.
- **Never** sort includes (`SortIncludes: false`) — it breaks OpenFOAM
  lnInclude symlink resolution order

### Naming conventions
| Kind | Convention | Example |
|---|---|---|
| Classes | `lowerCamelCase` | `solidBodyDisplacementLaplacianZoneFvMotionSolver` |
| Private members | `lowerCamelCase_` (trailing `_`) | `cellDisplacement_`, `diffusivityPtr_` |
| Functions / methods | `lowerCamelCase` | `curPoints()`, `updateMesh()` |
| Local variables | `lowerCamelCase` | `cellIDs`, `zoneID` |
| Constants / enums | `UPPER_CASE` | — |

### Types — always prefer OpenFOAM types
- `label` not `int` for indices and counts
- `scalar` not `double` for floating-point
- `word` not `std::string` for dictionary identifiers
- `vector`, `tensor`, `symmTensor` for geometric quantities
- `tmp<Field<T>>` for temporary field returns
- `autoPtr<T>` or `std::unique_ptr<T>` for sole-ownership pointers
- `Foam::List<T>` for MPI-serialisable lists; `std::vector<T>` is acceptable
  internally but convert via range constructor: `Foam::List<T>(v.begin(), v.end())`

### Runtime type registration (in `.C`, inside `namespace Foam`)
```cpp
namespace Foam
{
    defineTypeNameAndDebug(ClassName, 0);
    addToRunTimeSelectionTable(baseClass, ClassName, dictionary);
}
```

### Error handling
- `FatalErrorInFunction << "msg" << exit(FatalError);` — unrecoverable errors
- `FatalIOErrorInFunction(dict) << "msg" << exit(FatalIOError);` — I/O errors
- `WarningInFunction << "msg" << endl;` — non-fatal warnings
- `Info<< "msg" << endl;` — informational output
- Always check `findZoneID()`, `findObject<>()`, and similar lookup return values

### MPI / parallel
- Use `OPstream`/`IPstream` for explicit rank-to-rank sends
- `syncTools::syncPointList()` / `syncTools::syncFaceList()` for halo sync
- `returnReduce(val, sumOp<label>())` for global reductions
- Gather-scatter: populate `List<List<T>>[myProcNo()]`, then use
  `preciceAdapter::gather()` / `preciceAdapter::scatter()`

### Mathematical constants
```cpp
#include "mathematicalConstants.H"
Foam::constant::mathematical::pi   // NOT M_PI
```

## Shell scripts
- All `Allwmake`/`Allclean` use `#!/bin/sh` (POSIX sh, not bash)
- First line: `cd "${0%/*}" || exit`
- Clean scripts use `set -e -u`
- Lint with `shellcheck --exclude=SC1091`

## Git
- Do not commit build artifacts: `*.o`, `*.so`, `platforms/`, `lnInclude/`, `*.log`
- License for all new files: GPL-3.0-or-later
- The adapter follows `precice-openfoam-adapter/CONTRIBUTING.md` for PRs
  (per-PR changelog entries in `changelog-entries/`, pre-commit hooks required)
