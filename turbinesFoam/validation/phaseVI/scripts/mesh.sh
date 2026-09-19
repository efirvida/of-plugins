#!/bin/sh
# Generate and validate one Phase VI mesh.
#
# Usage: mesh.sh <coarse|fine|ultra> [--domain long|squat] [--speed 7] [--dir DIR]
#
# The case skeleton is rendered into runs/mesh-<mesh>[-<domain>]/ and meshed
# with blockMesh + topoSet + checkMesh. The checkMesh log is validated against
# the configured cell-count band, the 100 % hexahedral requirement and the
# Cartesian non-orthogonality tolerance.
#
# Exit codes:
#   0  mesh generated and accepted
#   2  unsupported argument
#   3  blockMesh/checkMesh failed or the mesh is outside the accepted band
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$here/.." && pwd)

usage() {
    sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
}

mesh=""
domain=long
speed=7
dir=""
while [ $# -gt 0 ]; do
    case "$1" in
        coarse|fine|ultra) mesh="$1" ;;
        --domain) domain="$2"; shift ;;
        --speed) speed="$2"; shift ;;
        --dir) dir="$2"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unsupported argument: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done

if [ -z "$mesh" ]; then
    echo "ERROR: mesh resolution required (coarse|fine|ultra)" >&2
    usage >&2
    exit 2
fi
case "$domain" in
    long|squat) ;;
    *) echo "ERROR: unsupported domain: $domain" >&2; exit 2 ;;
esac

if [ -z "$dir" ]; then
    dir="$root/runs/mesh-$mesh"
    if [ "$domain" = "squat" ]; then
        dir="$dir-squat"
    fi
fi

echo "Rendering $mesh case ($domain domain) into $dir"
python3 "$root/tools/generate_case.py" \
    --mesh "$mesh" --speed "$speed" --domain "$domain" --case-dir "$dir"

cd "$dir"

# Re-use an existing mesh (blockMesh is serial and D/48 takes minutes); the
# checkMesh validation below still runs on every invocation. topoSet is
# re-run whenever its dictionary is newer than the last topoSet log: the
# cellSet is part of the case contract (`selectionMode cellSet` in fvOptions)
# and a mesh built from an older topoSetDict silently lacks the current set.
if [ -d constant/polyMesh ] && grep -q "Mesh OK" log.checkMesh 2>/dev/null; then
    echo "Existing mesh found in $dir; skipping blockMesh and re-validating"
    if [ ! -f log.topoSet ] || [ system/topoSetDict -nt log.topoSet ]; then
        echo "topoSetDict changed since the last topoSet run; refreshing the cellSet"
        if ! topoSet > log.topoSet 2>&1; then
            cat log.topoSet >&2
            echo "ERROR: topoSet failed ($dir)" >&2
            exit 3
        fi
    fi
else
    if ! blockMesh > log.blockMesh 2>&1; then
        cat log.blockMesh >&2
        echo "ERROR: blockMesh failed ($dir)" >&2
        exit 3
    fi
    if ! topoSet > log.topoSet 2>&1; then
        cat log.topoSet >&2
        echo "ERROR: topoSet failed ($dir)" >&2
        exit 3
    fi
fi
if ! checkMesh > log.checkMesh 2>&1; then
    cat log.checkMesh >&2
    echo "ERROR: checkMesh failed ($dir)" >&2
    exit 3
fi

grep -q "Mesh OK" log.checkMesh || {
    echo "ERROR: checkMesh did not report 'Mesh OK' ($dir)" >&2
    exit 3
}

cells=$(awk '/cells:/ {print $2; exit}' log.checkMesh)
hexes=$(awk '/hexahedra:/ {print $2; exit}' log.checkMesh)
polyhedra=$(awk '/polyhedra:/ {print $2; exit}' log.checkMesh)
nonortho=$(awk '/Mesh non-orthogonality Max:/ {print $4; exit}' log.checkMesh)

if [ -z "$cells" ] || [ -z "$hexes" ] || [ -z "$polyhedra" ] || [ -z "$nonortho" ]; then
    echo "ERROR: could not read the mesh statistics from $dir/log.checkMesh" >&2
    exit 3
fi
if [ "$hexes" != "$cells" ] || [ "$polyhedra" != "0" ]; then
    echo "ERROR: mesh is not 100% hexahedral ($hexes/$cells hexes, $polyhedra polyhedra)" >&2
    exit 3
fi
if ! awk -v value="$nonortho" 'BEGIN { exit !(value <= 1e-10) }'; then
    echo "ERROR: max non-orthogonality $nonortho exceeds the Cartesian tolerance" >&2
    exit 3
fi

band=$(python3 "$root/tools/case_config.py" --target-cells "$mesh")
low=${band% *}
high=${band#* }
if [ "$cells" -lt "$low" ] || [ "$cells" -gt "$high" ]; then
    echo "ERROR: $cells cells is outside the configured band $low..$high" >&2
    exit 3
fi

echo "Mesh accepted: $cells cells, $hexes hexahedra, max non-orthogonality $nonortho"
