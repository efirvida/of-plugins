#!/usr/bin/env python
"""Tests for the standalone `nacelleSurfaceSource` case.

The case (`tests/nacelleSurface`) is rotor-less: a two-triangle plate with the
outward normal (1 0 0) sits in a uniform 0.1 m cubic mesh (10 x 6 x 6 cells), the
inlet is a uniform (1 0 0) velocity, dt = 0.1, cf is a constant 0.01 and
rho = 1. That makes the first (and only) force evaluation analytic:

    f_n   = h * (u . e_n) / dt            (Eq. 19, u^d = 0)
          = 0.1 * 1 / 0.1 = 1.0
    f_tau = 0.5 * cf * U^2                (Eq. 21)
          = 0.5 * 0.01 * 1 = 0.005
    F     = rho * (f_n + f_tau) * A       (area weight)
          = (1.0 + 0.005) * 0.04 = 0.0402 N on the body

The distributed `force.nacelle` field is the force on the flow, so its volume
integral is minus the node-force sum: -0.0402. The case runs the kernel over a
5 x 5 x 5 stencil (h = 0.1, support 2.5 h = 0.25 in each direction), which
straddles the x = 0.5 processor boundary of `decomposeParDict`, so the parallel
run exercises cross-rank interpolation and distribution.
"""

from __future__ import division, print_function

import glob
import os
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

CASE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "nacelleSurface"
)

# Analytic case constants (see tests/nacelleSurface/system)
H_CELL = 0.1
DELTA_T = 0.1
U_REF = 1.0
CF = 0.01
RHO = 1.0
N_NODES = 2
AREA_NODE = 0.02
AREA_TOTAL = N_NODES * AREA_NODE
REFERENCE_AREA = 0.04

F_N_NODE = H_CELL * U_REF / DELTA_T
F_TAU_NODE = 0.5 * CF * U_REF**2
F_TOTAL_X = RHO * (F_N_NODE + F_TAU_NODE) * AREA_TOTAL
F_N_TOTAL = RHO * F_N_NODE * AREA_TOTAL
F_TAU_TOTAL = RHO * F_TAU_NODE * AREA_TOTAL
CD = F_TOTAL_X / (0.5 * RHO * U_REF**2 * REFERENCE_AREA)

PERF_CSV = os.path.join("postProcessing", "nacelle", "nacelle.csv")
NODES_CSV = os.path.join("postProcessing", "nacelle", "nacelle_nodes.csv")

# Totals of the serial run, for the parallel comparison
SERIAL = {}


@pytest.fixture(scope="module")
def case():
    """Run each test from the case directory, restoring the cwd afterwards."""
    orig = os.getcwd()
    os.chdir(CASE_DIR)
    yield
    os.chdir(orig)


def run_allrun(*args):
    """Clean and run the case; print the solver log and fail on error."""
    subprocess.check_output("./Allclean")
    try:
        subprocess.check_output(["./Allrun"] + list(args))
    except subprocess.CalledProcessError:
        print("".join(open("log.pimpleFoam").readlines()[-200:]))
        raise AssertionError("Allrun failed: {}".format(" ".join(args)))


def read_total():
    """Last row of the total-force CSV as a Series."""
    df = pd.read_csv(PERF_CSV)
    return df.drop_duplicates("time", keep="last").iloc[-1]


def read_force_integral():
    """Volume integral of the distributed `force.nacelle` field (on the flow)."""
    dat = sorted(glob.glob("postProcessing/forceIntegral/*/volFieldValue.dat"))[-1]
    rows = [line for line in open(dat) if not line.startswith("#")]
    # volIntegrate writes "<time>\t(<x> <y> <z>)"
    vec = rows[-1].split()[-3:]
    return np.array([float(v.strip("()")) for v in vec])


def check_force_csv():
    """The total and per-node CSVs match the analytic expectations."""
    perf = read_total()

    # Total force on the body, Eq. 19 + Eq. 21 with the area weight (SI N)
    assert_allclose(perf.fx, F_TOTAL_X, rtol=1e-6, atol=1e-12)
    assert_allclose([perf.fy, perf.fz], [0.0, 0.0], atol=1e-12)
    assert_allclose(perf.f_n_mag, F_N_TOTAL, rtol=1e-6)
    assert_allclose(perf.f_tau_mag, F_TAU_TOTAL, rtol=1e-6)
    assert_allclose(perf.cd, CD, rtol=1e-6)

    # The split diagnostics add up to the total (both act along +x here)
    assert_allclose(perf.f_n_mag + perf.f_tau_mag, perf.fx, rtol=1e-6)

    # Per-node CSV: same node set, body frame, SI forces
    nodes = pd.read_csv(NODES_CSV)
    assert len(nodes) == N_NODES
    assert len(set(nodes.node)) == N_NODES
    assert_allclose(nodes.area, [AREA_NODE] * N_NODES, rtol=1e-12)
    assert_allclose(nodes.nx, [1.0] * N_NODES, rtol=1e-12)
    assert_allclose(nodes.ny, [0.0] * N_NODES, atol=1e-12)
    assert_allclose(nodes.nz, [0.0] * N_NODES, atol=1e-12)
    assert_allclose(
        nodes.fx, [RHO * (F_N_NODE + F_TAU_NODE) * AREA_NODE] * N_NODES,
        rtol=1e-6,
    )


def check_total_force_preserved():
    """Integral of the distributed field equals the sum of the node forces.

    The kernel is volume-normalised (partition of unity), so
    `sum_c forceField[c] * V_c == sum_i F_i`; the distributed field is the
    force on the flow (Eq. 18's negative sign), i.e. minus the body force.
    """
    integral = read_force_integral()
    assert_allclose(integral[0], -F_TOTAL_X / RHO, rtol=1e-6, atol=1e-14)
    assert_allclose(integral[1:], [0.0, 0.0], atol=1e-14)


def test_serial(case):
    """Standalone source runs rotor-less and writes the force CSVs."""
    run_allrun()

    log = open("log.pimpleFoam").read()
    assert "Selecting finite volume options type nacelleSurfaceSource" in log
    assert log.rstrip().endswith("End")

    check_force_csv()
    check_total_force_preserved()

    SERIAL["fx"] = read_total().fx
    SERIAL["integral"] = read_force_integral()[0]


def test_total_force_preserved(case):
    """Total-force preservation of the last completed run (serial or parallel)."""
    if not os.path.isfile(PERF_CSV):
        pytest.skip("no case output; run test_serial or test_parallel first")
    check_total_force_preserved()


def test_parallel(case):
    """Two ranks: node->cell resolution across the x = 0.5 processor boundary."""
    # The plate sits at x = 0.45 and the kernel support is 2.5 h = 0.25, so the
    # 5 x 5 x 5 stencil spans x in [0.2, 0.7] and crosses the split at x = 0.5
    assert os.path.isfile(os.path.join(CASE_DIR, "system", "decomposeParDict"))
    assert "numberOfSubdomains 2" in open(
        os.path.join(CASE_DIR, "system", "decomposeParDict")
    ).read()

    run_allrun("-parallel")

    log = open("log.pimpleFoam").read()
    assert "Finalising parallel run" in log
    assert os.path.isdir("processor0")
    assert os.path.isdir("processor1")

    check_force_csv()
    check_total_force_preserved()

    # Same totals as serial within tolerance (both are analytic here)
    if SERIAL:
        assert_allclose(read_total().fx, SERIAL["fx"], rtol=1e-9)
        assert_allclose(read_force_integral()[0], SERIAL["integral"], rtol=1e-9)


def _copy_case(dst):
    """Copy the case without generated run artifacts."""
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


def _run_with_geometry(case_dir, geometry):
    """Point the source at `geometry`, run the solver and return (code, log)."""
    fv_options = os.path.join(case_dir, "system", "fvOptions")
    text = open(fv_options).read()
    assert "geometry/plate.stl" in text
    with open(fv_options, "w") as fv:
        fv.write(text.replace("geometry/plate.stl", geometry))

    subprocess.check_call(["cp", "-rf", "0.org", "0"], cwd=case_dir)
    for cmd in ("blockMesh", "topoSet"):
        subprocess.check_call(
            cmd, cwd=case_dir, stdout=subprocess.DEVNULL
        )

    log_path = os.path.join(case_dir, "log.pimpleFoam")
    with open(log_path, "w") as log:
        proc = subprocess.run(
            ["pimpleFoam"],
            cwd=case_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    return proc.returncode, open(log_path).read()


def test_missing_stl_aborts(tmp_path):
    """A missing STL raises FatalError instead of running silently."""
    case_dir = _copy_case(str(tmp_path / "missing"))
    code, log = _run_with_geometry(case_dir, "geometry/doesNotExist.stl")
    assert code != 0
    assert "FOAM FATAL ERROR" in log
    assert "not found" in log


def test_empty_stl_aborts(tmp_path):
    """An empty STL raises FatalError instead of running silently."""
    case_dir = _copy_case(str(tmp_path / "empty"))
    code, log = _run_with_geometry(case_dir, "geometry/empty.stl")
    assert code != 0
    assert "FOAM FATAL ERROR" in log
    assert "contains no triangles" in log


def test_corrupt_stl_aborts(tmp_path):
    """An unparseable STL raises FatalError instead of running silently."""
    case_dir = _copy_case(str(tmp_path / "corrupt"))
    with open(os.path.join(case_dir, "geometry", "corrupt.stl"), "w") as stl:
        stl.write("this is not an STL file\n")
    code, log = _run_with_geometry(case_dir, "geometry/corrupt.stl")
    assert code != 0
    assert "FOAM FATAL ERROR" in log
