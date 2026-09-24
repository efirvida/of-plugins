#!/usr/bin/env python
"""Tests for the Du-Selig rotational augmentation (W1).

Two layers:

``test_du_selig_reference_values`` is a pure-Python re-implementation of the
Du-Selig Eqs. 9-12 that reproduces the hand-computed U13 sample and the near-tip
``fL -> -1/(2*pi)`` limit, and pins the C++ constants/exponents and the rendered
configuration so the reference cannot silently drift.

The remaining tests drive ``tests/rotationalAugmentation/rotor``, a minimal
one-blade ``axialFlowTurbineALSource`` fixture with an S809 polar. Its
``elementData`` endpoints are r = 0.05 m and r = 0.45 m of a 0.45 m rotor, so
the two element centres land at r = 0.15 m and r = 0.35 m (chord 0.04 m, no
twist, TSR 3). ``system/fvOptions`` is the default-off baseline; the ``.on``, ``.off``,
``.asm``, ``.on-fast``, ``.nogeom`` and ``.badmodel`` variants install through
``./Allrun`` on a throwaway copy of the case, so the committed files are never
mutated.

The oracle recomputes the expected corrected coefficients from the CSV's own
angle of attack using a Python replica of ``profileData``'s interpolation (the
fixture omits the ``Re`` key, so the static polar is the raw committed table).
"""

from __future__ import division, print_function

import os
import re
import shutil
import subprocess

import glob

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

CASE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "rotationalAugmentation", "rotor"
)
# One CSV per element, under the run's start-time directory. A restarted run
# writes under its restart time, so the lookup globs every time directory.
ELEMENT_CSV_GLOB = os.path.join(
    "postProcessing", "actuatorLineElements", "*", "turbine.blade1.element{}.csv"
)
POLAR_DAT = os.path.join(CASE_DIR, "system", "S809.dat")

# Fixture constants (tests/rotationalAugmentation/rotor/system/fvOptions)
U_REF = 10.0
ROTOR_RADIUS = 0.45
ROOT_RADIUS = 0.05
CHORD = 0.04
N_ELEMENTS = 2
TSR_ON = 3.0
TSR_PRE_STALL = 50.0

# Element centres from the fixture's `elementData` endpoints (0.05 m and
# 0.45 m). With one geometry segment split into two elements, `createElements`
# places them at the quarter and three-quarter points: r = 0.15 m and 0.35 m.
ELEMENT_RADIUS = {0: 0.15, 1: 0.35}


def _read_polar():
    """Parse the committed S809 polar table into (alpha, cl, cd)."""
    alpha, cl, cd = [], [], []
    for line in open(POLAR_DAT):
        line = line.strip()
        if not line.startswith("("):
            continue
        vals = [float(v) for v in line.strip("()").split()]
        alpha.append(vals[0])
        cl.append(vals[1])
        cd.append(vals[2])
    return np.array(alpha), np.array(cl), np.array(cd)


ALPHA_TABLE, CL_TABLE, CD_TABLE = _read_polar()


def _interpolate(x_new, x_old, y_old):
    """Replicate ``profileData::interpolate`` (nearest index + linear span)."""
    x_old = np.asarray(x_old, dtype=float)
    y_old = np.asarray(y_old, dtype=float)
    index = int(np.argmin(np.abs(x_new - x_old)))
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
    """Replicate ``profileData::subList`` over the angle-of-attack list."""
    return [full[i] for i in range(len(ALPHA_TABLE))
            if start <= ALPHA_TABLE[i] <= stop]


def _zero_lift_constants():
    """``zeroLiftAngleOfAttack`` (deg) and ``zeroLiftDragCoeff`` of the polar."""
    cl = _sublist(-10, 10, CL_TABLE)
    alpha = _sublist(-10, 10, ALPHA_TABLE)
    cd = _sublist(-10, 10, CD_TABLE)
    return _interpolate(0.0, cl, alpha), _interpolate(0.0, cl, cd)


ALPHA0_DEG, CD0 = _zero_lift_constants()


def static_cl(alpha_deg):
    return _interpolate(alpha_deg, ALPHA_TABLE, CL_TABLE)


def static_cd(alpha_deg):
    return _interpolate(alpha_deg, ALPHA_TABLE, CD_TABLE)


def du_selig_factors(c_over_r, r_over_r, lam, a=1.0, b=1.0, d=1.0):
    """Du-Selig Eqs. 11-12 factors ``fL`` and ``fD``."""
    x_l = (d / lam) * r_over_r
    x_d = (d / (2.0 * lam)) * r_over_r
    q_l = c_over_r ** x_l
    q_d = c_over_r ** x_d
    f_l = (1.0 / (2.0 * np.pi)) * (
        (1.6 * c_over_r ** a - q_l) / (0.1267 * b + q_l) - 1.0
    )
    f_d = (1.0 / (2.0 * np.pi)) * (
        (1.6 * c_over_r ** a - q_d) / (0.1267 * b + q_d) - 1.0
    )
    return f_l, f_d


def du_selig_corrected(f_l, f_d, cl_p, cl_2d, cd_2d):
    """Du-Selig Eqs. 9-10 given the factors and the potential-flow lift."""
    return cl_2d + f_l * (cl_p - cl_2d), cd_2d - f_d * (cd_2d - CD0)


def du_selig(c_over_r, r_over_r, lam, alpha_deg, cl_2d, cd_2d,
             a=1.0, b=1.0, d=1.0):
    """Du-Selig Eqs. 9-12 with ``CL,p = 2*pi*(alpha - alpha0)``."""
    f_l, f_d = du_selig_factors(c_over_r, r_over_r, lam, a, b, d)
    alpha = np.deg2rad(alpha_deg)
    alpha_0 = np.deg2rad(ALPHA0_DEG)
    cl_p = 2.0 * np.pi * (alpha - alpha_0)
    cl_3d, cd_3d = du_selig_corrected(f_l, f_d, cl_p, cl_2d, cd_2d)
    return f_l, f_d, cl_p, cl_3d, cd_3d


def _expected_corrected(alpha_deg, radius, tsr):
    """Expected (cl, cd) for the fixture operating point at ``alpha_deg``."""
    omega = tsr * U_REF / ROTOR_RADIUS
    omega_r = omega * ROTOR_RADIUS
    lam = omega_r / np.sqrt(U_REF ** 2 + omega_r ** 2)
    cl_2d = static_cl(alpha_deg)
    cd_2d = static_cd(alpha_deg)
    _, _, _, cl_3d, cd_3d = du_selig(
        CHORD / radius, ROTOR_RADIUS / radius, lam, alpha_deg, cl_2d, cd_2d
    )
    return cl_3d, cd_3d, cl_2d, cd_2d


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


def _run_case(case_dir, *args):
    """Run ``./Allrun`` in ``case_dir``; fail with the log tail on error."""
    proc = subprocess.run(
        ["./Allrun"] + list(args),
        cwd=case_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    if proc.returncode != 0:
        log_path = os.path.join(case_dir, "log.pimpleFoam")
        tail = open(log_path).readlines()[-200:] if os.path.isfile(log_path) else []
        pytest.fail(
            "Allrun failed ({}) exit {}:\n{}".format(
                " ".join(args), proc.returncode, "".join(tail) or proc.stdout
            )
        )


def _run_case_expect_failure(case_dir, *args):
    """Run ``./Allrun`` and return (returncode, log text) without failing."""
    proc = subprocess.run(
        ["./Allrun"] + list(args),
        cwd=case_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    log_path = os.path.join(case_dir, "log.pimpleFoam")
    log = open(log_path).read() if os.path.isfile(log_path) else proc.stdout
    return proc.returncode, log


def _read_element(case_dir, element=0):
    """Read one element's CSV across every written time, sorted by time."""
    paths = sorted(glob.glob(os.path.join(
        case_dir, ELEMENT_CSV_GLOB.format(element)
    )))
    assert paths, "no element CSV for element {} under {}".format(element, case_dir)
    df = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
    return df.sort_values("time").reset_index(drop=True)


def _initial_row(case_dir, element=0):
    """The first written row for an element.

    ``pimpleFoam`` increments the clock before the first force evaluation, so
    the element CSV starts at the first time step, not at ``t = 0``; the first
    row is the one evaluated on the uniform initial field and is therefore the
    deterministic operating point the reference reproduces.
    """
    df = _read_element(case_dir, element)
    return df.iloc[0]


def _read_log(case_dir):
    return open(os.path.join(case_dir, "log.pimpleFoam")).read()


def _set_control(case_dir, start_from=None, end_time=None):
    path = os.path.join(case_dir, "system", "controlDict")
    text = open(path).read()
    if start_from is not None:
        text = re.sub(r"startFrom\s+\w+;", "startFrom       %s;" % start_from, text)
    if end_time is not None:
        text = re.sub(r"endTime\s+[\d.]+;", "endTime         %s;" % end_time, text)
    open(path, "w").write(text)


@pytest.fixture(scope="module")
def default_case(tmp_path_factory):
    case_dir = str(tmp_path_factory.mktemp("ra-default") / "case")
    _copy_case(case_dir)
    _run_case(case_dir)
    return case_dir


@pytest.fixture(scope="module")
def off_case(tmp_path_factory):
    case_dir = str(tmp_path_factory.mktemp("ra-off") / "case")
    _copy_case(case_dir)
    _run_case(case_dir, "-off")
    return case_dir


@pytest.fixture(scope="module")
def on_case(tmp_path_factory):
    case_dir = str(tmp_path_factory.mktemp("ra-on") / "case")
    _copy_case(case_dir)
    _run_case(case_dir, "-on")
    return case_dir


@pytest.fixture(scope="module")
def asm_case(tmp_path_factory):
    case_dir = str(tmp_path_factory.mktemp("ra-asm") / "case")
    _copy_case(case_dir)
    _run_case(case_dir, "-asm")
    return case_dir


@pytest.fixture(scope="module")
def prestall_case(tmp_path_factory):
    case_dir = str(tmp_path_factory.mktemp("ra-prestall") / "case")
    _copy_case(case_dir)
    _run_case(case_dir, "-on-fast")
    return case_dir


@pytest.fixture(scope="module")
def nogeom_case(tmp_path_factory):
    case_dir = str(tmp_path_factory.mktemp("ra-nogeom") / "case")
    _copy_case(case_dir)
    _run_case(case_dir, "-nogeom")
    return case_dir


@pytest.fixture(scope="module")
def cylinder_case(tmp_path_factory):
    case_dir = str(tmp_path_factory.mktemp("ra-cylinder") / "case")
    _copy_case(case_dir)
    _run_case(case_dir, "-cylinder")
    return case_dir


# --------------------------------------------------------------------------- #
# Pure-Python reference
# --------------------------------------------------------------------------- #

def test_du_selig_reference_values():
    """Eqs. 9-12 reproduce the hand sample, the near-tip limit and the C++."""
    # Design/spec hand sample: c/r ~ 0.271, R/r ~ 2.27, Lambda ~ 0.95,
    # CL,p ~ 3.1 and the diagnosed static CL,2D in 0.598..0.73 -> 1.11..1.21.
    f_l, _ = du_selig_factors(0.271, 2.27, 0.95)
    assert f_l == pytest.approx(0.2036, rel=2e-3)
    cl_3d, _ = du_selig_corrected(f_l, 0.0, 3.1, 0.73, 0.0)
    assert cl_3d == pytest.approx(1.213, rel=1e-2)
    cl_3d_control, _ = du_selig_corrected(f_l, 0.0, 3.1, 0.598, 0.0)
    assert 1.10 <= cl_3d_control <= 1.22

    # Near-tip limit: as c/r -> 0, fL -> -1/(2*pi) (recorded, never clamped).
    f_l_tip, _ = du_selig_factors(1.0e-4, 20.0, 0.95)
    assert f_l_tip == pytest.approx(-1.0 / (2.0 * np.pi), abs=0.01)
    assert f_l_tip < 0.0

    # Exponent form is (d/Lambda)(R/r), never d*Lambda*R/r: the latter would
    # make fL larger for a larger Lambda, the former smaller. Pin the direction
    # so a text-dump-style exponent cannot silently pass.
    f_l_low_lambda, _ = du_selig_factors(0.271, 2.27, 0.5)
    assert f_l_low_lambda > f_l  # smaller Lambda -> larger exponent -> larger fL


def test_cpp_constants_and_config_pinned():
    """The Python reference cannot drift from the C++ and the rendered dict."""
    src = open(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "src", "fvOptions", "actuatorLineSource", "actuatorLineElement",
        "actuatorLineElement.C",
    )).read()
    # Equation constants and the verified exponent forms.
    assert "1.6*Foam::pow(cOverR, a_)" in src
    assert "0.1267*b_" in src
    assert "(d_/lambda)*ROverR" in src
    assert "(d_/(2.0*lambda))*ROverR" in src
    assert "liftCoefficient_ += fL*(CLp - liftCoefficient_)" in src
    assert "dragCoefficient_ -= fD*(dragCoefficient_ - CD0)" in src
    # No invented clamp of fL.
    assert "max(fL" not in src and "min(fL" not in src

    rendered = open(os.path.join(
        CASE_DIR, "system", "fvOptions.on"
    )).read()
    assert "rotationalAugmentation" in rendered
    assert "active          on;" in rendered
    assert "model           DuSelig;" in rendered
    for key in ("a", "b", "d"):
        assert re.search(r"^\s+%s\s+1;" % key, rendered, re.M)


def test_order_relative_to_dynamic_stall():
    """The hook runs after the static lookup and before dynamic stall.

    The combined Du-Selig + Leishman-Beddoes model is not claimed validated
    (design §2.5), so the ordering contract is pinned structurally: the
    augmentation call must sit between ``lookupCoefficients()`` and the
    ``dynamicStall_->correct`` call in ``calculateForce``. W2 exercises the LB
    path at runtime; the augmentation-first order is fixed here.
    """
    src = open(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "src", "fvOptions", "actuatorLineSource", "actuatorLineElement",
        "actuatorLineElement.C",
    )).read()
    i_lookup = src.index("lookupCoefficients();")
    i_hook = src.index("if (rotationalAugmentationActive_)")
    i_augment = src.index("correctRotationalAugmentation();")
    i_dynamic_stall = src.index("dynamicStall_->correct")
    assert i_lookup < i_hook < i_augment < i_dynamic_stall


# --------------------------------------------------------------------------- #
# Integration
# --------------------------------------------------------------------------- #

def test_default_off_byte_identical(default_case, off_case):
    """No block and explicit ``active off`` produce identical element output."""
    for element in range(N_ELEMENTS):
        default_paths = sorted(glob.glob(os.path.join(
            default_case, ELEMENT_CSV_GLOB.format(element)
        )))
        off_paths = sorted(glob.glob(os.path.join(
            off_case, ELEMENT_CSV_GLOB.format(element)
        )))
        assert default_paths and off_paths
        assert open(default_paths[0], "rb").read() == open(off_paths[0], "rb").read()

    # Non-circular gate: the default path is exactly the static polar, i.e. the
    # new code never runs when the block is absent.
    for element in range(N_ELEMENTS):
        row = _initial_row(default_case, element)
        assert row["cl"] == pytest.approx(static_cl(row["alpha_deg"]), rel=1e-4)
        assert row["cd"] == pytest.approx(static_cd(row["alpha_deg"]), rel=1e-4)


def test_correction_applied(default_case, on_case):
    """The active block rewrites cl/cd exactly per Eqs. 9-12."""
    for element in range(N_ELEMENTS):
        radius = ELEMENT_RADIUS[element]
        row = _initial_row(on_case, element)
        alpha = row["alpha_deg"]
        cl_3d, cd_3d, cl_2d, _ = _expected_corrected(alpha, radius, TSR_ON)
        assert row["cl"] == pytest.approx(cl_3d, rel=2e-4)
        assert row["cd"] == pytest.approx(cd_3d, rel=2e-4)
        # The correction changed the output from the static polar.
        assert row["cl"] != pytest.approx(cl_2d, rel=1e-6)
        # Inboard deep stall gains lift (fL > 0); outboard loses it (fL < 0).
        if element == 0:
            assert row["cl"] > cl_2d
        else:
            assert row["cl"] < cl_2d
        assert row["cl"] != pytest.approx(
            _initial_row(default_case, element)["cl"]
        )


def test_radial_geometry_keys(default_case, on_case):
    """``createElements`` injects ``radius``/``rotorRadius`` additively.

    The injected station must equal the comparison tool's identity
    ``r = rootRadius + rootDistance*(rotorRadius - rootRadius)``; inverting the
    CSV's ``root_dist`` recovers the element centre, and the corrected ``cl``
    matches the reference only when that radius is the one the element used.
    """
    for element in range(N_ELEMENTS):
        row = _initial_row(on_case, element)
        injected = ROOT_RADIUS + row["root_dist"] * (ROTOR_RADIUS - ROOT_RADIUS)
        assert injected == pytest.approx(ELEMENT_RADIUS[element], rel=1e-9)
        cl_3d, _, _, _ = _expected_corrected(row["alpha_deg"], injected, TSR_ON)
        assert row["cl"] == pytest.approx(cl_3d, rel=2e-4)

    # No block -> no radial keys are injected and the static polar is untouched.
    for element in range(N_ELEMENTS):
        row = _initial_row(default_case, element)
        assert row["cl"] == pytest.approx(static_cl(row["alpha_deg"]), rel=1e-4)
        assert row["cd"] == pytest.approx(static_cd(row["alpha_deg"]), rel=1e-4)


def test_pre_stall_unchanged(prestall_case):
    """At a pre-stall angle the lift correction vanishes (spec scenario)."""
    # Outboard element (r = 0.35 m) at TSR 50 sits at alpha ~ 1.47 deg, where
    # CL,2D ~ CL,p.
    row = _initial_row(prestall_case, 1)
    alpha = row["alpha_deg"]
    assert alpha < 5.0
    cl_3d, _, cl_2d, _ = _expected_corrected(alpha, ELEMENT_RADIUS[1], TSR_PRE_STALL)
    assert cl_3d == pytest.approx(cl_2d, abs=1e-3)
    assert row["cl"] == pytest.approx(cl_3d, rel=2e-4)


def test_absent_keys_keep_old_dict(nogeom_case):
    """Active block without geometry warns and skips; no abort."""
    for element in range(N_ELEMENTS):
        row = _initial_row(nogeom_case, element)
        assert row["cl"] == pytest.approx(static_cl(row["alpha_deg"]), rel=1e-4)
        assert row["cd"] == pytest.approx(static_cd(row["alpha_deg"]), rel=1e-4)
    log = _read_log(nogeom_case)
    assert "rotationalAugmentation active but radius/rotorRadius absent" in log


def test_degenerate_profile_skipped(cylinder_case):
    """A two-point `cylinder` root profile is skipped, never a crash.

    The Phase VI case mixes a degenerate root `cylinder` (only +/-180 deg) with
    the S809 lifting sections. The augmentation's zero-lift lookup interpolates
    over [-10, 10] deg, which is empty for the cylinder; without the guard that
    indexes an empty list and segfaults. The lifting element must still be
    corrected.
    """
    # Element 0 (cylinder): static table only, cl = 0 and cd = 1.1, no crash.
    row0 = _initial_row(cylinder_case, 0)
    assert row0["cl"] == pytest.approx(0.0, abs=1e-9)
    assert row0["cd"] == pytest.approx(1.1, rel=1e-6)

    # Element 1 (S809): the correction is applied exactly per Eqs. 9-12.
    row1 = _initial_row(cylinder_case, 1)
    cl_3d, cd_3d, cl_2d, _ = _expected_corrected(
        row1["alpha_deg"], ELEMENT_RADIUS[1], TSR_ON
    )
    assert row1["cl"] == pytest.approx(cl_3d, rel=2e-4)
    assert row1["cd"] == pytest.approx(cd_3d, rel=2e-4)
    assert row1["cl"] != pytest.approx(cl_2d, rel=1e-6)

    # The guard is structural too: the element skips on the profile query.
    src = open(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "src", "fvOptions", "actuatorLineSource", "actuatorLineElement",
        "actuatorLineElement.C",
    )).read()
    assert "profileData_.hasZeroLiftReference()" in src
    profile_src = open(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "src", "fvOptions", "actuatorLineSource", "actuatorLineElement",
        "profileData", "profileData.C",
    )).read()
    assert "bool Foam::profileData::hasZeroLiftReference()" in profile_src


def test_unknown_model_rejected(tmp_path):
    """An unregistered model fails loudly (FatalIOError)."""
    case_dir = str(tmp_path / "case")
    _copy_case(case_dir)
    returncode, log = _run_case_expect_failure(case_dir, "-badmodel")
    assert returncode != 0
    assert "Unknown rotationalAugmentation model" in log
    assert "notRegistered" in log


def test_all_models_inherit(on_case, asm_case):
    """ALM and the no-mesh surface element share the same corrected chain."""
    for case_dir in (on_case, asm_case):
        for element in range(N_ELEMENTS):
            row = _initial_row(case_dir, element)
            alpha = row["alpha_deg"]
            cl_3d, cd_3d, cl_2d, _ = _expected_corrected(
                alpha, ELEMENT_RADIUS[element], TSR_ON
            )
            assert row["cl"] == pytest.approx(cl_3d, rel=2e-4)
            assert row["cd"] == pytest.approx(cd_3d, rel=2e-4)
            assert row["cl"] != pytest.approx(cl_2d, rel=1e-6)


def test_parallel_matches_serial(on_case, tmp_path):
    """``mpirun -np 2`` reproduces the serial corrected rows."""
    case_dir = str(tmp_path / "case")
    _copy_case(case_dir)
    _run_case(case_dir, "-on", "-parallel")
    for element in range(N_ELEMENTS):
        serial = _initial_row(on_case, element)
        parallel = _initial_row(case_dir, element)
        assert_allclose(
            [parallel["cl"], parallel["cd"], parallel["alpha_deg"],
             parallel["root_dist"]],
            [serial["cl"], serial["cd"], serial["alpha_deg"],
             serial["root_dist"]],
            rtol=1e-6,
            atol=1e-9,
        )


def test_restart_reproduces_scalars(tmp_path):
    """A restart from a written time reproduces the corrected scalars."""
    # Continuous run to t = 0.003.
    continuous = str(tmp_path / "continuous" / "case")
    _copy_case(continuous)
    _set_control(continuous, end_time="0.003")
    _run_case(continuous, "-on")
    continuous_df = _read_element(continuous, 0)
    assert continuous_df["time"].iloc[-1] == pytest.approx(0.003)
    continuous_row = continuous_df.iloc[-1]

    # Run to t = 0.002, then restart from latestTime to t = 0.003.
    restarted = str(tmp_path / "restarted" / "case")
    _copy_case(restarted)
    _set_control(restarted, end_time="0.002")
    _run_case(restarted, "-on")
    _set_control(restarted, start_from="latestTime", end_time="0.003")
    proc = subprocess.run(
        ["pimpleFoam"],
        cwd=restarted,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    if proc.returncode != 0:
        pytest.fail("restart pimpleFoam failed:\n" + proc.stdout[-2000:])
    restart_df = _read_element(restarted, 0)
    assert restart_df["time"].iloc[-1] == pytest.approx(0.003)
    restart_row = restart_df.iloc[-1]

    assert restart_row["alpha_deg"] == pytest.approx(
        continuous_row["alpha_deg"], rel=1e-4
    )
    assert restart_row["cl"] == pytest.approx(continuous_row["cl"], rel=1e-3)
    assert restart_row["cd"] == pytest.approx(continuous_row["cd"], rel=1e-3)
