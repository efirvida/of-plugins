#!/bin/sh
# Verify the environment and the generated case for the Phase VI validation.
#
# Accepts OpenFOAM v2506 (repository toolchain) and v2412.
#
# Exit codes:
#   0  environment OK and generated case is current
#   3  environment problem (OpenFOAM version, missing tool or library)
#   4  generated case is missing or stale
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$here/.." && pwd)

version="${WM_PROJECT_VERSION:-}"
case "$version" in
    v2506|v2412) ;;
    *)
        echo "ERROR: load OpenFOAM v2506 or v2412 first (found '${version:-unset}')." >&2
        echo "  module load openfoam/v2506_openmpi-4.1.4_gnu" >&2
        exit 3
        ;;
esac

missing=0
for command in blockMesh topoSet checkMesh pimpleFoam decomposePar mpirun python3; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $command" >&2
        missing=1
    fi
done
if [ "$missing" -ne 0 ]; then
    exit 3
fi

library="${FOAM_USER_LIBBIN:-}/libturbinesFoam.so"
if [ ! -f "$library" ]; then
    echo "ERROR: libturbinesFoam.so not found in FOAM_USER_LIBBIN." >&2
    echo "  Build it with: cd turbinesFoam && ./Allwmake" >&2
    exit 3
fi

if ! python3 "$root/tools/generate_case.py" --check; then
    echo "ERROR: generated case is stale; re-render with tools/generate_case.py" >&2
    exit 4
fi

echo "Environment OK: OpenFOAM $version"
echo "turbinesFoam library: $library"
echo "Generated case: consistent with config/case.yaml"
