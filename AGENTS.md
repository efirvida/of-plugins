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

The adapter is an FSI-only fork of the upstream preCICE adapter: coupling
modules live in `precice-openfoam-adapter/modules/FSI/` and
`modules/generic/` (`Stress` and `DisplacementDelta` modules were removed).

Built library names (loaded via `libs (...)` in `controlDict`):
`libsolidBodyDisplacementLaplacianZoneFvMotionSolver`,
`libdynamicOversetZoneDisplacementFvMesh`, `libfsiOmega`,
`libpreciceAdapterFunctionObject`. Root `README.md` documents case setup
and the FSI motion model.

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

Requires a loaded OpenFOAM environment (plugin `Allclean` scripts refuse to
run without `WM_PROJECT_VERSION`). The adapter's `Allwmake` uses
`pkg-config libprecice` automatically and emits warnings (not errors) when
the `.pc` file is missing. Environment variables:
- `PRECICE_OPENFOAM_CFLAGS` — extra preprocessor flags (e.g. `-DADAPTER_DEBUG_MODE`)
- `PRECICE_OPENFOAM_TARGET_DIR` — install destination (default `$FOAM_USER_LIBBIN`)

## Testing

**No automated test suite.** Validation is manual — run a real OpenFOAM case
that loads the compiled library. When adding new code:
- `./Allwmake` must exit 0; the adapter's script additionally fails on any
  `error:` line in `wmake.log` and runs `ldd -r`, checking for undefined symbols
- For MPI-related changes, test in parallel with `mpirun -np N`

Adapter system regression tests run on an external system, only when a
maintainer adds the `trigger-system-tests` label to the PR (see
`precice-openfoam-adapter/CONTRIBUTING.md`).

## CI/CD (precice-openfoam-adapter only)

Workflows live in `precice-openfoam-adapter/.github/workflows/`:
- `build.yml` — builds against OpenFOAM **v2512** (preCICE nightly container) on push/PR
- `pre-commit.yml` — runs pre-commit hooks on push to main/develop and all PRs
- `check-format.yml` — clang-format v11 on all `.C`/`.H` files
- `check-shell.yml` — shellcheck on all shell scripts
- `build-custom.yml` — manual multi-version build (via `build-custom.input` + `act`)
- `system-tests.yaml` — label-triggered regression tests

## Code formatting

The authoritative style is `precice-openfoam-adapter/.clang-format` (clang-format 11).
```sh
cd precice-openfoam-adapter
sh tools/format-code.sh               # format all .C/.H files
pip install pre-commit && pre-commit install
pre-commit run -va                    # apply hooks retrospectively
```

Pre-commit hooks (`.pre-commit-config.yaml`): clang-format v11, shellcheck
(`--exclude=SC1091`), markdownlint (`docs/*.md` only), preCICE config checks.

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

### OpenFOAM v2506+ gotchas (CI builds v2512)
- `Foam::List<T>` **cannot** be constructed from STL iterators (unsupported
  constructors) — use `setSize()` + element copy loops (see `Interface.C`)
- `dictionary::lookup()` is deprecated — use `get<word>()`
- `Pstream::scatterList`/`gatherList` are deprecated and broken in v2506 —
  use `preciceAdapter::gather()`/`scatter()` from `Utilities.H`, or
  `OPstream`/`IPstream`

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
- `Foam::List<T>` for MPI-serialisable lists; `std::vector<T>` acceptable
  internally (convert manually — see gotchas above)

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
- `syncTools::syncPointList()` / `syncTools::syncFaceList()` for halo sync
- `returnReduce(val, sumOp<label>())` for global reductions
- Gather-scatter: populate `List<List<T>>[myProcNo()]`, then use
  `preciceAdapter::gather()` / `preciceAdapter::scatter()`

### Mathematical constants
```cpp
#include "mathematicalConstants.H"
Foam::constant::mathematical::pi   // NOT M_PI
```

## Documentation conventions
- Root `CHANGELOG.md` documents each notable change in
  `Files:` / `Problem:` / `Fix:` format — update it alongside code changes,
  and document new solver options in root `README.md`
- Adapter changes: add `changelog-entries/<PR-number>.md` (never edit
  `precice-openfoam-adapter/CHANGELOG.md` directly — entries are merged at release)

## Shell scripts
- All `Allwmake`/`Allclean` use `#!/bin/sh` (POSIX sh, not bash)
- First line: `cd "${0%/*}" || exit`
- Clean scripts use `set -e -u`
- Lint with `shellcheck --exclude=SC1091`

## Git
- Do not commit build artifacts: `*.o`, `*.so`, `platforms/`, `lnInclude/`, `*.log`
- License for all new files: GPL-3.0-or-later
- Adapter PRs follow `precice-openfoam-adapter/CONTRIBUTING.md`
  (changelog entries, pre-commit hooks, system tests via label)
