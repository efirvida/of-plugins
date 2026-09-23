#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Stage the committed Phase VI blade STL into a prepared run directory.

`runPhaseVI.sh -m asm-mesh` renders the run directory from `config/case.yaml`,
installs the `fvOptions.ASM-MESH` twin and then calls this helper. The helper
verifies the committed STL against the `stl_sha256` recorded in the committed
geometry metadata, copies it to `constant/triSurface/phaseVI_blade.stl` (the
case-relative path the ASM-mesh twin references as `surfaceGeometry`) and
prints the staged sha256 on stdout, so the runner can record it in `run.json`
and the comparison can assert the fair-comparison input hash.

Staging is what makes a surface-model run self-contained: restarts and
relocated runs resolve `surfaceGeometry` inside their own directory, and a
mismatch aborts preparation instead of running with a stale surface.

Usage:
    stage_blade_stl.py --run-dir DIR [--stl PATH] [--metadata PATH]

Exit codes:
    0  the STL was staged and its sha256 matches the metadata
    3  staging failure (missing source/metadata or sha256 mismatch); the
       caller aborts preparation
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STL = ROOT.parents[1] / "geometry" / "stl" / "phaseVI_blade.stl"
DEFAULT_METADATA = ROOT.parents[1] / "geometry" / "metadata" / "phaseVI_blade.json"
# Case-relative destination: what the ASM-mesh twin references.
STAGED_RELATIVE = Path("constant") / "triSurface" / "phaseVI_blade.stl"
EXIT_STAGING = 3


class StageError(RuntimeError):
    """A staging failure that aborts run preparation."""


def sha256_file(path: Path) -> str:
    """Hex sha256 of a file, read in chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expected_stl_sha256(metadata: Path) -> str:
    """The committed STL sha256 recorded in the geometry metadata."""
    try:
        payload = json.loads(metadata.read_text(encoding="utf-8"))
    except OSError as exc:
        raise StageError(f"cannot read blade metadata {metadata}: {exc}") from exc
    except ValueError as exc:
        raise StageError(
            f"blade metadata {metadata} is not valid JSON: {exc}"
        ) from exc
    digest = payload.get("stl_sha256") if isinstance(payload, dict) else None
    if not isinstance(digest, str) or not digest:
        raise StageError(f"blade metadata {metadata} does not record stl_sha256")
    return digest


def stage_blade_stl(
    stl: Path = DEFAULT_STL,
    metadata: Path = DEFAULT_METADATA,
    run_dir: Path | None = None,
) -> str:
    """Stage `stl` into `<run_dir>/constant/triSurface/phaseVI_blade.stl`.

    The source is checked against the metadata `stl_sha256` before the copy
    and the staged copy is checked again afterwards; the staged sha256 is
    returned. Raises `StageError` on any missing input or hash mismatch.
    """
    if run_dir is None:
        raise StageError("a run directory is required")
    stl = Path(stl)
    metadata = Path(metadata)
    if not stl.is_file():
        raise StageError(f"committed blade STL not found: {stl}")
    expected = expected_stl_sha256(metadata)
    actual = sha256_file(stl)
    if actual != expected:
        raise StageError(
            f"blade STL sha256 mismatch: {stl} is {actual}, "
            f"metadata records {expected}"
        )
    destination = Path(run_dir) / STAGED_RELATIVE
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(stl, destination)
    staged = sha256_file(destination)
    if staged != expected:
        raise StageError(
            f"staged STL sha256 mismatch: {destination} is {staged}, "
            f"metadata records {expected}"
        )
    return staged


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="prepared run directory that receives constant/triSurface/",
    )
    parser.add_argument(
        "--stl",
        type=Path,
        default=DEFAULT_STL,
        help="committed blade STL (default: geometry/stl/phaseVI_blade.stl)",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=DEFAULT_METADATA,
        help="committed blade metadata carrying stl_sha256",
    )
    args = parser.parse_args(argv)
    try:
        staged = stage_blade_stl(args.stl, args.metadata, args.run_dir)
    except StageError as exc:
        print(f"staging error: {exc}", file=sys.stderr)
        return EXIT_STAGING
    # Bare hash on stdout so the runner can capture it for run.json; the
    # human-readable location goes to stderr.
    print(staged)
    print(
        f"staged {Path(args.run_dir) / STAGED_RELATIVE} (sha256 {staged})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
