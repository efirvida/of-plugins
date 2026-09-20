# SPDX-License-Identifier: GPL-3.0-or-later
"""Test-suite configuration.

The upstream turbinesFoam suites drive OpenFOAM tutorial cases (`./Allrun`)
and require a loaded OpenFOAM environment plus the built `libturbinesFoam.so`.
The NREL Phase VI suites are pure Python. When OpenFOAM is not loaded, skip
the solver-driven modules so that

    python -m pytest turbinesFoam/tests -q

is meaningful on a login node without an OpenFOAM installation.
"""

import os
import shutil

SOLVER_DRIVEN = (
    "test_aftal.py",
    "test_aftal_asm.py",
    "test_al.py",
    "test_asm.py",
    "test_cftal.py",
    "test_libs.py",
    "test_nacelle.py",
)


def _openfoam_available():
    if not os.environ.get("WM_PROJECT_VERSION"):
        return False
    return any(
        shutil.which(name) for name in ("blockMesh", "simpleFoam", "pimpleFoam")
    )


collect_ignore = [] if _openfoam_available() else list(SOLVER_DRIVEN)
