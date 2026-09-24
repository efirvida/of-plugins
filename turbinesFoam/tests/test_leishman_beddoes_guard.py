#!/usr/bin/env python
"""Tests for the Leishman-Beddoes K1/K2 singularity guard (W2).

Two layers:

``test_well_conditioned_matches_solve`` and ``test_singular_matrix_falls_back``
are pure-Python references for the analytic 2x2 (Cramer) solve and the
determinant guard in ``LeishmanBeddoes::calcK1K2``. They replicate the normal
equations from the fixture's raw S809 polar and pin the C++ expressions so the
reference cannot silently drift. ``test_blast_radius_one_method`` and
``test_no_fit_retuning`` pin the two structural spec scenarios.

``test_guard_run_passes_crash_point`` drives ``tests/leishmanBeddoes/rotor``, a
minimal one-blade ``axialFlowTurbineALSource`` fixture. Element 0 is a cylinder
polar with no station in the K1/K2 fit range [0.5, 25] deg, so its normal-
equation matrix is identically zero and the guard must take the documented
fallback; element 1 is the S809 polar with a well-conditioned matrix, so the
analytic branch is exercised in the same run. The fixture supplies explicit
``S1List``/``S2List`` so the S1/S2 fit (whose ranges are subsets of the empty
K1/K2 interval) is skipped and ``calcK1K2`` is the only Leishman-Beddoes fit
exercised. The former crash point is the first PIMPLE iteration at
``t = 0.008 s`` (``deltaT = 0.008``).
"""

from __future__ import division, print_function

import glob
import math
import os
import re
import shutil
import subprocess

import numpy as np
import pytest

CASE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "leishmanBeddoes", "rotor"
)
S809_DAT = os.path.join(CASE_DIR, "system", "S809.dat")
FV_OPTIONS = os.path.join(CASE_DIR, "system", "fvOptions")

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_TURBINESFOAM = os.path.dirname(_TESTS_DIR)
_LEISHMANBEDDOES_C = os.path.join(
    _TURBINESFOAM, "src", "fvOptions", "actuatorLineSource",
    "actuatorLineElement", "dynamicStallModels", "LeishmanBeddoes",
    "LeishmanBeddoes.C",
)
_PROFILEDATA_C = os.path.join(
    _TURBINESFOAM, "src", "fvOptions", "actuatorLineSource",
    "actuatorLineElement", "profileData", "profileData.C",
)
_DYNAMIC_STALL_DIR = os.path.join(
    _TURBINESFOAM, "src", "fvOptions", "actuatorLineSource",
    "actuatorLineElement", "dynamicStallModels",
)

# LeishmanBeddoes constructor defaults (LeishmanBeddoes.C:550,553) and the
# OpenFOAM scalar VSMALL used by the guard.
K0 = 1.0e-6
CM_FIT_EXPONENT = 2
VSMALL = 1.0e-300
SMALL = 1.0e-15

FALLBACK_WARNING = "singular/ill-conditioned Leishman-Beddoes K1/K2 fit matrix"
FALLBACK_LINE = "falling back to K1 = K2 = 0"
FIT_PATH_MARKER = "path = well-conditioned)"


# --------------------------------------------------------------------------- #
# profileData / LeishmanBeddoes pure-Python replica
# --------------------------------------------------------------------------- #

def _read_s809():
    """Parse the committed S809 polar into (alpha_deg, cl, cd, cm)."""
    alpha, cl, cd, cm = [], [], [], []
    for line in open(S809_DAT):
        line = line.strip()
        if not line.startswith("("):
            continue
        vals = [float(v) for v in line.strip("()").split()]
        alpha.append(vals[0])
        cl.append(vals[1])
        cd.append(vals[2])
        cm.append(vals[3])
    return alpha, cl, cd, cm


ALPHA_TABLE, CL_TABLE, CD_TABLE, CM_TABLE = _read_s809()


def _interpolate(x_new, x_old, y_old):
    """Replicate ``profileData::interpolate`` (nearest index + linear span)."""
    x_old = list(x_old)
    y_old = list(y_old)
    index = min(range(len(x_old)), key=lambda i: abs(x_new - x_old[i]))
    if x_new < x_old[index]:
        if index == 0:
            index_p, index_m = 1, 0
        else:
            index_p, index_m = index, index - 1
    elif x_new > x_old[index]:
        if index == len(x_old) - 1:
            index_p, index_m = len(x_old) - 1, len(x_old) - 2
        else:
            index_p, index_m = index + 1, index
    else:
        return float(y_old[index])
    return float(
        y_old[index_m]
        + (y_old[index_p] - y_old[index_m])
        / (x_old[index_p] - x_old[index_m])
        * (x_new - x_old[index_m])
    )


def _sublist(start, stop, full):
    """Replicate ``profileData::subList`` (inclusive on the raw table order)."""
    return [full[i] for i in range(len(ALPHA_TABLE))
            if start <= ALPHA_TABLE[i] <= stop]


def _cn(alpha_deg, cl, cd):
    """``profileData::convertToCN``."""
    r = math.radians(alpha_deg)
    return cl * math.cos(r) + cd * math.sin(r)


def _static_stall_angle():
    """Replicate ``profileData::calcStaticStallAngle`` (first cd-slope break)."""
    for a in ALPHA_TABLE:
        if 2 < a < 30 and (
            _interpolate(a + 1.0, ALPHA_TABLE, CD_TABLE)
            - _interpolate(a, ALPHA_TABLE, CD_TABLE)
        ) > 0.03:
            return a
    raise AssertionError("S809 polar has no static stall break")


def _normal_coeff_slope(alpha_high):
    """Replicate ``profileData::calcNormalCoeffSlope`` (slope of cn vs alpha)."""
    a_list = _sublist(0, alpha_high, ALPHA_TABLE)
    cn_list = [
        _cn(ALPHA_TABLE[i], CL_TABLE[i], CD_TABLE[i])
        for i in range(len(ALPHA_TABLE)) if 0 <= ALPHA_TABLE[i] <= alpha_high
    ]
    n = len(a_list)
    x = [math.radians(v) for v in a_list]
    sx = sum(x)
    sy = sum(cn_list)
    sxx = sum(v * v for v in x)
    sxy = sum(x[i] * cn_list[i] for i in range(n))
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def _zero_lift_angle():
    cl_10 = _sublist(-10, 10, CL_TABLE)
    alpha_10 = _sublist(-10, 10, ALPHA_TABLE)
    return _interpolate(0.0, cl_10, alpha_10)


def _cn_to_f(cn_list, alpha_rad_list, cn_alpha, alpha0_rad, limit=False):
    """Replicate ``LeishmanBeddoes::cnToF`` (``calcK1K2`` uses ``limit=False``)."""
    f_list = []
    for i in range(len(cn_list)):
        alpha_rad = alpha_rad_list[i]
        if alpha_rad == 0.0:
            alpha_rad = SMALL
        f = (
            2.0 * math.sqrt(abs(cn_list[i] / (cn_alpha * abs(alpha_rad - alpha0_rad))))
            - 1.0
        ) ** 2
        if limit:
            if f >= 1:
                f = 1.0 - SMALL
            if f <= 0:
                f = SMALL
        f_list.append(f)
    return f_list


def _calc_k1k2_normal_equations():
    """The exact ``calcK1K2`` normal-equation matrix for the raw S809 polar.

    Mirrors ``calcStaticStallAngle``/``calcNormalCoeffSlope``/``cnToF`` and the
    ``[0.5, 25]`` deg fit range. The fixture's ``Re 1e6`` makes the runtime
    lists Reynolds-corrected; this reference is the raw-polar representative
    point, which is what the pure-Python oracle pins.
    """
    alpha_high = 0.5 * _static_stall_angle()
    cn_alpha = _normal_coeff_slope(alpha_high)
    alpha0_rad = math.radians(_zero_lift_angle())

    alpha_deg = _sublist(0.5, 25, ALPHA_TABLE)
    cn_list = _sublist(
        0.5, 25,
        [_cn(ALPHA_TABLE[i], CL_TABLE[i], CD_TABLE[i])
         for i in range(len(ALPHA_TABLE))],
    )
    cm_list = _sublist(0.5, 25, CM_TABLE)
    alpha_rad = [math.radians(v) for v in alpha_deg]
    f = _cn_to_f(cn_list, alpha_rad, cn_alpha, alpha0_rad)

    m = CM_FIT_EXPONENT
    sin_pi_f = [math.sin(math.pi * fi ** m) for fi in f]
    A = np.array([
        [sum((1.0 - fi) ** 2 for fi in f), sum(sin_pi_f)],
        [sum(sin_pi_f[i] * (1.0 - f[i]) for i in range(len(f))),
         sum(s ** 2 for s in sin_pi_f)],
    ])
    b = np.array([
        sum(cm_list[i] / cn_list[i] * (1.0 - f[i]) for i in range(len(f)))
        - K0 * sum(1.0 - fi for fi in f),
        sum(cm_list[i] / cn_list[i] * sin_pi_f[i] for i in range(len(f)))
        - K0 * sum(sin_pi_f),
    ])
    return A, b


def _guard(A, b):
    """Replicate the C++ determinant guard and return (K1, K2, det, path)."""
    det = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    scale = max(abs(A[0, 0] * A[1, 1]), abs(A[0, 1] * A[1, 0])) + VSMALL
    if abs(det) > 1.0e-12 * scale:
        k1 = (b[0] * A[1, 1] - A[0, 1] * b[1]) / det
        k2 = (A[0, 0] * b[1] - b[0] * A[1, 0]) / det
        return k1, k2, det, "well-conditioned"
    return 0.0, 0.0, det, "fallback"


def _read_source(path):
    return open(path).read()


# --------------------------------------------------------------------------- #
# Pure-Python guard oracle
# --------------------------------------------------------------------------- #

def test_well_conditioned_matches_solve():
    """The analytic solve equals the unguarded (LU) solve to round-off."""
    A, b = _calc_k1k2_normal_equations()
    assert A.shape == (2, 2)

    k1, k2, det, path = _guard(A, b)
    scale = max(abs(A[0, 0] * A[1, 1]), abs(A[0, 1] * A[1, 0])) + VSMALL
    # The representative fit is genuinely well-conditioned, so the guard must
    # take the analytic branch (never the fallback).
    assert path == "well-conditioned"
    assert abs(det) > 1.0e-12 * scale
    assert abs(det) / scale > 0.1

    # The unguarded reference: numpy's LU solve of the same normal equations.
    reference = np.linalg.solve(A, b)
    assert k1 == pytest.approx(reference[0], rel=1e-12, abs=1e-14)
    assert k2 == pytest.approx(reference[1], rel=1e-12, abs=1e-14)
    # Pin the representative values so a silent fixture/replica drift fails.
    assert k1 == pytest.approx(-0.126573, rel=1e-5)
    assert k2 == pytest.approx(0.077436, rel=1e-5)

    # The C++ must carry exactly these analytic expressions.
    src = _read_source(_LEISHMANBEDDOES_C)
    assert "const scalar det = A[0][0]*A[1][1] - A[0][1]*A[1][0];" in src
    assert (
        "Foam::max(mag(A[0][0]*A[1][1]), mag(A[0][1]*A[1][0])) + VSMALL" in src
    )
    assert "if (mag(det) > 1.0e-12*scale)" in src
    assert (
        "K1_ = (A.source()[0]*A[1][1] - A[0][1]*A.source()[1])/det;" in src
    )
    assert (
        "K2_ = (A[0][0]*A.source()[1] - A.source()[0]*A[1][0])/det;" in src
    )


def test_singular_matrix_falls_back():
    """A singular/dependent matrix yields K1 = K2 = 0 with no abort."""
    # The cylinder fixture's empty [0.5, 25] deg fit range gives A = b = 0.
    A_zero = np.zeros((2, 2))
    b_zero = np.zeros(2)
    k1, k2, det, path = _guard(A_zero, b_zero)
    assert (k1, k2) == (0.0, 0.0)
    assert det == 0.0
    assert path == "fallback"

    # A rank-deficient (row-proportional) matrix is also detected.
    A_dep = np.array([[1.0, 2.0], [2.0, 4.0]])
    b_dep = np.array([1.0, 2.0])
    k1, k2, det, path = _guard(A_dep, b_dep)
    assert (k1, k2) == (0.0, 0.0)
    assert det == 0.0
    assert path == "fallback"

    # The fallback is documented in the C++ and records the determinant.
    src = _read_source(_LEISHMANBEDDOES_C)
    assert FALLBACK_WARNING in src
    assert FALLBACK_LINE in src
    assert "WarningInFunction" in src
    assert "K1_ = 0.0;" in src
    assert "K2_ = 0.0;" in src


def test_blast_radius_one_method():
    """Only ``LeishmanBeddoes::calcK1K2`` is modified; others are untouched."""
    sources = glob.glob(os.path.join(_DYNAMIC_STALL_DIR, "**", "*.C"), recursive=True)
    assert sources
    guarded = [p for p in sources if FALLBACK_LINE in _read_source(p)]
    assert guarded == [_LEISHMANBEDDOES_C]

    # The unguarded solve is gone from the LB model...
    lb_src = _read_source(_LEISHMANBEDDOES_C)
    assert "A.solve()" not in lb_src
    assert "List<scalar> sol = A.solve();" not in lb_src

    # ...and ``profileData::calcNormalCoeffSlope`` still uses it unchanged.
    assert "normalCoeffSlope_ = A.solve()[1];" in _read_source(_PROFILEDATA_C)


def test_no_fit_retuning():
    """The guard does not alter the fit inputs (spec: no fit retuning)."""
    src = _read_source(_LEISHMANBEDDOES_C)
    # Fit inputs and the documented fallback value are unchanged.
    assert 'cmFitExponent_(coeffs_.lookupOrDefault("cmFitExponent", 2))' in src
    assert "K0_(1e-6)" in src
    # The fallback removes the moment contribution; it does not invent a fit.
    assert "K1_ = 0.0;\n        K2_ = 0.0;" in src


def test_fixture_cylinder_fit_range_is_empty():
    """Pin the integration premise: the cylinder has no station in [0.5, 25]."""
    text = _read_source(FV_OPTIONS)
    cylinder = re.search(
        r"cylinder\s*\{(.*?)\n            \}", text, re.S
    )
    assert cylinder, "cylinder profileData block not found"
    alphas = [
        float(m.group(1))
        for m in re.finditer(r"\(\s*(-?\d+(?:\.\d+)?)\s", cylinder.group(1))
    ]
    assert alphas, "no cylinder polar rows parsed"
    assert not [a for a in alphas if 0.5 <= a <= 25.0]


# --------------------------------------------------------------------------- #
# Case harness
# --------------------------------------------------------------------------- #

def _copy_case(dst):
    """Copy the fixture without generated run artifacts."""
    shutil.copytree(CASE_DIR, dst)
    for name in os.listdir(dst):
        path = os.path.join(dst, name)
        generated = (
            name in ("0", "postProcessing")
            or name.startswith("processor")
            or name.startswith("log.")
            or (name[0].isdigit() and name != "0.org")
        )
        if not generated:
            continue
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
    return dst


def _run_case(case_dir):
    """Run ``./Allrun``; fail with the log tail on error."""
    proc = subprocess.run(
        ["./Allrun"],
        cwd=case_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    if proc.returncode != 0:
        log_path = os.path.join(case_dir, "log.pimpleFoam")
        tail = open(log_path).readlines()[-200:] if os.path.isfile(log_path) else []
        pytest.fail(
            "Allrun failed exit {}:\n{}".format(
                proc.returncode, "".join(tail) or proc.stdout
            )
        )


def _read_log(case_dir):
    return open(os.path.join(case_dir, "log.pimpleFoam")).read()


def _enable_leishman_beddoes_debug(case_dir):
    """Turn on the ``LeishmanBeddoes`` DebugSwitch for the throwaway copy.

    OpenFOAM reads ``DebugSwitches`` from the case ``system/controlDict``
    (``TimeIO.C``), which is what makes the well-conditioned ``if (debug)``
    record observable in the run log.
    """
    path = os.path.join(case_dir, "system", "controlDict")
    text = open(path).read()
    assert "DebugSwitches" not in text
    text += "\nDebugSwitches\n{\n    LeishmanBeddoes 1;\n}\n"
    open(path, "w").write(text)


@pytest.fixture(scope="module")
def guard_case(tmp_path_factory):
    case_dir = str(tmp_path_factory.mktemp("lb-guard") / "case")
    _copy_case(case_dir)
    _run_case(case_dir)
    return case_dir


@pytest.fixture(scope="module")
def guard_debug_case(tmp_path_factory):
    case_dir = str(tmp_path_factory.mktemp("lb-guard-debug") / "case")
    _copy_case(case_dir)
    _enable_leishman_beddoes_debug(case_dir)
    _run_case(case_dir)
    return case_dir


# --------------------------------------------------------------------------- #
# Integration
# --------------------------------------------------------------------------- #

def test_guard_run_passes_crash_point(guard_case):
    """The run passes t = 0.008 s with no ``Singular Matrix`` abort."""
    log = _read_log(guard_case)
    assert "FOAM FATAL ERROR" not in log
    assert "Singular Matrix" not in log
    # The former crash point and the following step are both reached.
    assert "Time = 0.008" in log
    assert "Time = 0.016" in log
    # The fallback path is recorded with its determinant (det = 0 for the
    # cylinder's all-zero normal-equation matrix).
    assert FALLBACK_WARNING in log
    assert FALLBACK_LINE in log
    assert re.search(r"det = 0\b", log)

    # Both elements were evaluated (element CSV written).
    csvs = glob.glob(os.path.join(
        guard_case, "postProcessing", "actuatorLineElements", "*",
        "turbine.blade1.element*.csv",
    ))
    assert len(csvs) >= 2


def test_well_conditioned_path_recorded(guard_debug_case):
    """With debug on, the well-conditioned analytic path is recorded."""
    log = _read_log(guard_debug_case)
    assert "FOAM FATAL ERROR" not in log
    assert "Singular Matrix" not in log
    assert "Leishman-Beddoes K1/K2 analytic solve (det = " in log
    assert FIT_PATH_MARKER in log
