#!/usr/bin/env python
"""Integration tests for the blade surface distributor (W1.7).

Fixtures
--------

``tests/bladeSurface`` is a rotor-less standalone ``actuatorLineSource`` case:
a 0.1 m cubic mesh (26 x 26 x 26 cells), uniform inflow ``(1 0 0)``, one
cylinder element whose chord line lies in the two-triangle plate STL
``geometry/plate.stl``, and a ``forceIntegral`` functionObject on
``force.blade``. ``system/fvOptions`` is the delivered default (no surface);
``system/fvOptions.surface`` and ``system/fvOptions.surface.gaussian``
activate the distributor (cosine and Gaussian kernel). Every test runs on a
throwaway copy of the case and installs the variant through ``./Allrun``, so
the committed case files are never mutated.

``tests/bladeSurfaceAFTAL`` is the minimal rotating fixture used by the
rotation/moment tests: one blade, two geometry rows (one element), no
hub/tower, ``writeNodePerf true`` and ``rho 1.0``, with a multi-triangle plate
placed 0.09 m outboard of the element station. The spanwise offset is
deliberate: a purely chordwise offset leaves the rotor-axis moment unchanged
(the section force has no spanwise component), so only the spanwise offset
lets the moment test tell the D8 surface moment from the element-position
moment. The surface CSVs are written with 12 significant digits, because the
moment oracle (``rtol 1e-6`` against ``turbine.csv``) cannot be represented
with the delivered six-significant-digit element/turbine output.

Analytic oracle at the first step (t = 0.1, uniform field, cylinder profile
``cd = 1.1``, element velocity zero)::

    area = chord * span = 0.2 * 0.2 = 0.04 m2
    F    = 0.5 * area * cd * |U|^2 = 0.022 N on the blade (+x)

``force.blade`` is the force on the flow, so its volume integral is minus the
blade force. The surface distribution is a partition of unity
(``sum_cells w = 1``) and its integral equals the node/patch total to machine
precision. The delivered strip projection is a truncated Gaussian whose
discrete sum is 0.99978 of the element force in this fixture (the support is
not clipped by the domain), so the default-path test uses a documented 1e-3
tolerance instead of an exact equality.
"""

from __future__ import division, print_function

import glob
import os
import re
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

CASE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "bladeSurface"
)
AFTAL_CASE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "bladeSurfaceAFTAL"
)

STATION_CSV = os.path.join("postProcessing", "bladeSurface", "blade.surface.csv")
NODE_CSV = os.path.join(
    "postProcessing", "bladeSurface", "blade.surface_nodes.csv"
)
DISTRIBUTION_CSV = os.path.join(
    "postProcessing", "bladeSurface", "blade.surface_distribution.csv"
)
ELEMENT_CSV = os.path.join(
    "postProcessing", "actuatorLineElements", "0", "blade.element0.csv"
)
AFTAL_STATION_CSV = os.path.join(
    "postProcessing", "bladeSurface", "turbine.blade1.surface.csv"
)
AFTAL_NODE_CSV = os.path.join(
    "postProcessing", "bladeSurface", "turbine.blade1.surface_nodes.csv"
)
AFTAL_ELEMENT_CSV = os.path.join(
    "postProcessing", "actuatorLineElements", "0", "turbine.blade1.element0.csv"
)
TURBINE_CSV = os.path.join("postProcessing", "turbines", "0", "turbine.csv")

STATION_COLUMNS = [
    "time",
    "station",
    "root_dist",
    "area",
    "force_x",
    "force_y",
    "force_z",
    "c_ref_n",
    "c_ref_t",
    "f_ref_n",
    "f_ref_t",
]
NODE_COLUMNS = [
    "time",
    "node",
    "x",
    "y",
    "z",
    "nx",
    "ny",
    "nz",
    "fx",
    "fy",
    "fz",
    "area",
    "station",
    "chord_fraction",
]
DISTRIBUTION_COLUMNS = [
    "time",
    "nodes",
    "candidates",
    "mean_candidates",
    "max_candidates",
    "seconds",
]

# The per-addSup instrumentation line (D6/W4.1). The format is pinned here so
# the spec's "Measurement output present" scenario is an executable oracle.
DISTRIBUTION_RE = re.compile(
    r"Blade surface distribution 'blade\.surface': "
    r"nodes (\d+), candidates (\d+), mean ([0-9eE.+-]+), "
    r"max (\d+), seconds ([0-9eE.+-]+)"
)

# Analytic fixture constants (see tests/bladeSurface/system)
N_TIMES = 2
U_REF = 1.0
CD = 1.1
CHORD = 0.2
SPAN = 0.2
ELEMENT_AREA = CHORD * SPAN
F_ELEMENT = 0.5 * ELEMENT_AREA * CD * U_REF**2
PLATE_AREA = 0.04
N_NODES = 2
NODE_AREA = PLATE_AREA / N_NODES
# The triSurface face areas carry the STL vertex rounding (~2e-7 relative);
# the nominal plate area is the exact 0.04 m2 of the analytic oracle.
AREA_RTOL = 1e-6
ELEMENT_STATION = 0.5
ELEMENT_ROOT_DIST = 0.5

# Triangle centroids in the canonical generation frame (face order of the STL)
NODE_STATION = {0: 0.5333333333333333, 1: 0.4666666666666667}
NODE_CHORD_FRACTION = {0: 0.6666666666666666, 1: 0.3333333333333333}

# Fixture mesh size (tests/bladeSurface/system/blockMeshDict): 26^3 uniform
# 0.1 m cells. The naive distribution would scan every local cell per node
# (N_NODES * N_LOCAL_CELLS); the bounded query visits only the kernel support.
N_LOCAL_CELLS = 26**3

# Discrete sum of the delivered element strip projection in this fixture
# (truncated Gaussian, unclipped support): the default-path integral is this
# fraction of the element total.
STRIP_KERNEL_SUM = 0.999783

# Discrete normalization of the Gaussian ablation in this fixture (D7): the
# truncated Gaussian's discrete volume sum, measured as the ratio of the field
# integral to the patch total.
GAUSSIAN_KERNEL_SUM = 0.9995295

# AFTAL fixture (tests/bladeSurfaceAFTAL): one blade, one cylinder element at
# station 0.35 m, rho 1.0, per-node CSV enabled. The turbine runs about +x with
# freeStreamVelocity 10 and tipSpeedRatio 6 (rotorRadius 0.45), and injects the
# construction frame origin (0 0 0), span +z, chord -y (the element azimuthal
# direction), so the surface placement basis is fixed by the fixture.
AFTAL_N_NODES = 8
ROTOR_RADIUS = 0.45
U_INF = 10.0
TIP_SPEED_RATIO = 6.0
OMEGA = TIP_SPEED_RATIO * U_INF / ROTOR_RADIUS
FRONTAL_AREA = np.pi * ROTOR_RADIUS**2
TORQUE_DENOMINATOR = 0.5 * FRONTAL_AREA * ROTOR_RADIUS * U_INF**2
ROTOR_AXIS = np.array([1.0, 0.0, 0.0])


def _copy_case(src, dst):
    """Copy a case without generated run artifacts."""
    shutil.copytree(src, dst)
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
    """Run `./Allrun` in `case_dir`; fail with the log tail on error."""
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


def _read_log(case_dir):
    return open(os.path.join(case_dir, "log.pimpleFoam")).read()


def _read_force_integral(case_dir):
    """[(time, vector)] rows of the force.blade volume integral (on the flow)."""
    dat = sorted(
        glob.glob(
            os.path.join(case_dir, "postProcessing", "forceIntegral", "*",
                         "volFieldValue.dat")
        )
    )[-1]
    rows = []
    for line in open(dat):
        if line.startswith("#"):
            continue
        parts = line.split()
        vec = np.array([float(v.strip("()")) for v in parts[1:4]])
        rows.append((float(parts[0]), vec))
    return rows


def _read_station(case_dir):
    return pd.read_csv(os.path.join(case_dir, STATION_CSV))


def _read_nodes(case_dir):
    return pd.read_csv(os.path.join(case_dir, NODE_CSV))


def _read_distribution(case_dir):
    return pd.read_csv(os.path.join(case_dir, DISTRIBUTION_CSV))


def _read_element(case_dir):
    return pd.read_csv(os.path.join(case_dir, ELEMENT_CSV))


def _read_aftal_nodes(case_dir):
    return pd.read_csv(os.path.join(case_dir, AFTAL_NODE_CSV))


def _read_aftal_element(case_dir):
    return pd.read_csv(os.path.join(case_dir, AFTAL_ELEMENT_CSV))


def _read_turbine(case_dir):
    return pd.read_csv(os.path.join(case_dir, TURBINE_CSV))


def _body_basis(span):
    """Reproduce ``surfaceSamplerBase::createBodyFrame`` for a body axis.

    The completion is the C++ one: e2 = e1 ^ z, falling back to e1 ^ y when
    that is parallel, and e3 = e1 ^ e2. The columns are the body axes in the
    global frame, so ``global = basis @ body`` (the fixtures' bodyOrigin equals
    their surfaceOrigin, the origin).
    """
    e1 = np.asarray(span, dtype=float)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(e1, [0.0, 0.0, 1.0])
    if np.linalg.norm(e2) < 1e-12:
        e2 = np.cross(e1, [0.0, 1.0, 0.0])
    e2 = e2 / np.linalg.norm(e2)
    e3 = np.cross(e1, e2)
    return np.column_stack([e1, e2, e3])


def _to_global(values, basis):
    """Map body-frame rows to the global frame (fixture origin at (0 0 0))."""
    return (basis @ np.asarray(values, dtype=float).T).T


def _axis_rotation(radians):
    """Rotation matrix about the +x rotor axis (the ``rotateVector`` formula)."""
    c, s = np.cos(radians), np.sin(radians)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def _row_vector(row, columns):
    return row[columns].to_numpy(dtype=float)


@pytest.fixture(scope="module")
def surface_case(tmp_path_factory):
    """A completed `-surface` run of the standalone fixture."""
    case_dir = str(tmp_path_factory.mktemp("bladeSurface-surface") / "case")
    _copy_case(CASE_DIR, case_dir)
    _run_case(case_dir, "-surface")
    return case_dir


@pytest.fixture(scope="module")
def default_case(tmp_path_factory):
    """A completed default run (no surfaceGeometry) of the fixture."""
    case_dir = str(tmp_path_factory.mktemp("bladeSurface-default") / "case")
    _copy_case(CASE_DIR, case_dir)
    _run_case(case_dir)
    return case_dir


@pytest.fixture(scope="module")
def gaussian_case(tmp_path_factory):
    """A completed `-gaussian` run of the standalone fixture (D7 ablation)."""
    case_dir = str(tmp_path_factory.mktemp("bladeSurface-gaussian") / "case")
    _copy_case(CASE_DIR, case_dir)
    _run_case(case_dir, "-gaussian")
    return case_dir


@pytest.fixture(scope="module")
def aftal_case(tmp_path_factory):
    """A completed run of the minimal rotating AFTAL fixture."""
    case_dir = str(tmp_path_factory.mktemp("bladeSurfaceAFTAL") / "case")
    _copy_case(AFTAL_CASE_DIR, case_dir)
    _run_case(case_dir)
    return case_dir


def test_serial_surface(surface_case):
    """Partition of unity: the field sees each element force exactly once."""
    log = _read_log(surface_case)
    assert "Blade surface sampler: 2 nodes over 1 element patches" in log
    assert "kernel cosine" in log
    assert "Blade surface source 'blade.surface': 2 nodes, 1 element patches" in log

    station = _read_station(surface_case)
    assert len(station) == N_TIMES

    # First step, uniform inflow: exactly the analytic cylinder drag
    first = station.iloc[0]
    assert_allclose(
        _row_vector(first, ["force_x", "force_y", "force_z"]),
        [F_ELEMENT, 0.0, 0.0],
        rtol=1e-12,
        atol=1e-15,
    )

    # The distributed field is the force on the flow (minus the blade force)
    # and the kernel is a partition of unity: integral == -patch total. The
    # tolerance is the 6-significant-digit rounding of the written CSV.
    integral = _read_force_integral(surface_case)
    assert len(integral) == N_TIMES
    for time, vec in integral:
        patch = _row_vector(
            station[station.time == time].iloc[-1],
            ["force_x", "force_y", "force_z"],
        )
        assert_allclose(vec, -patch, rtol=1e-5, atol=1e-12)

    # The element CSV stays the BEM value: the patch total equals it exactly
    # (suppressed strip projection, not added to it)
    element = _read_element(surface_case)
    for time, _ in integral:
        element_total = _row_vector(
            element[element.time == time].iloc[-1], ["fx", "fy", "fz"]
        )
        patch = _row_vector(
            station[station.time == time].iloc[-1],
            ["force_x", "force_y", "force_z"],
        )
        assert_allclose(patch, element_total, rtol=1e-5, atol=1e-12)


def test_default_path_projects_strips(default_case):
    """No surface: the strip projection still reaches the field."""
    log = _read_log(default_case)
    assert "Blade surface sampler" not in log
    assert "Blade surface source" not in log
    assert not os.path.isfile(os.path.join(default_case, STATION_CSV))

    element = _read_element(default_case)
    total = _row_vector(element.iloc[0], ["fx", "fy", "fz"])
    assert_allclose(total, [F_ELEMENT, 0.0, 0.0], rtol=1e-12, atol=1e-15)

    # projectElementForce was never injected: the key only reaches an element
    # when a surface is configured or the user sets it, and a suppressed
    # element would leave the field integral at zero. The integral is the
    # element total to within the delivered strip projection's discrete
    # truncation (STRIP_KERNEL_SUM).
    time, vec = _read_force_integral(default_case)[0]
    assert time == 0.1
    assert np.linalg.norm(vec) > 0.0
    assert_allclose(-vec, total, rtol=1e-3, atol=1e-12)
    assert_allclose(-vec[0] / total[0], STRIP_KERNEL_SUM, rtol=1e-3)


def test_suppressed_element_still_reports(default_case, surface_case):
    """The suppressed element computes and reports; only the field changes."""
    default = _read_element(default_case)
    surface = _read_element(surface_case)

    for frame in (default, surface):
        assert len(frame) == N_TIMES
        assert np.all(
            np.isfinite(
                frame[["fx", "fy", "fz", "cl", "cd"]].to_numpy(dtype=float)
            )
        )
    assert np.linalg.norm(_row_vector(surface.iloc[-1], ["fx", "fy", "fz"])) > 0

    # The first PIMPLE solve of the first step sees the same uniform field in
    # both runs, so the first element row is identical
    assert list(default.columns) == list(surface.columns)
    assert_array_equal(
        surface.iloc[0].to_numpy(dtype=float),
        default.iloc[0].to_numpy(dtype=float),
    )

    # Later rows are finite and physical; they differ because the surface
    # distribution changed the inflow (the strip projection did not)
    assert not np.allclose(
        surface.iloc[-1].to_numpy(dtype=float),
        default.iloc[-1].to_numpy(dtype=float),
        rtol=1e-3,
        atol=1e-12,
    )

    # The field integral is the surface total, not surface + strips
    time, vec = _read_force_integral(surface_case)[-1]
    station = _read_station(surface_case)
    patch = _row_vector(
        station[station.time == time].iloc[-1],
        ["force_x", "force_y", "force_z"],
    )
    assert_allclose(vec, -patch, rtol=1e-5, atol=1e-12)


def test_station_and_node_csv(surface_case):
    """Both CSV schemas parse and carry the documented node metadata."""
    station = _read_station(surface_case)
    assert list(station.columns) == STATION_COLUMNS
    assert len(station) == N_TIMES
    assert station.time.nunique() == N_TIMES
    assert np.all(np.isfinite(station.to_numpy(dtype=float)))

    # One element, so one row per time step at the element station and the
    # element root-distance convention
    assert_allclose(station.station, ELEMENT_STATION, rtol=0, atol=1e-12)
    assert_allclose(station.root_dist, ELEMENT_ROOT_DIST, rtol=0, atol=1e-12)
    assert_allclose(station.area, PLATE_AREA, rtol=AREA_RTOL)
    assert np.all(station.area > 0)

    nodes = _read_nodes(surface_case)
    assert list(nodes.columns) == NODE_COLUMNS
    assert len(nodes) == N_TIMES * N_NODES
    assert np.all(np.isfinite(nodes.to_numpy(dtype=float)))

    for time in station.time:
        per_time = nodes[nodes.time == time]
        assert sorted(per_time.node) == [0, 1]
        assert_allclose(per_time.area, NODE_AREA, rtol=AREA_RTOL)

        # Chord fraction is LE-based and inside [0, 1]
        assert np.all(per_time.chord_fraction >= 0.0)
        assert np.all(per_time.chord_fraction <= 1.0)

        # Body frame contract: the body x-axis is the construction span, so
        # the node station equals its body x coordinate
        assert_allclose(per_time.station, per_time.x, rtol=1e-12, atol=1e-12)

        for _, row in per_time.iterrows():
            assert_allclose(
                row.station, NODE_STATION[row.node], rtol=1e-5, atol=1e-12
            )
            assert_allclose(
                row.chord_fraction,
                NODE_CHORD_FRACTION[row.node],
                rtol=1e-5,
                atol=1e-12,
            )

    # Node stations are monotone in radius (strictly increasing, STL order)
    unique_stations = np.sort(nodes.station.unique())
    assert_allclose(
        unique_stations, sorted(NODE_STATION.values()), rtol=1e-5, atol=1e-12
    )
    assert np.all(np.diff(unique_stations) > 0.0)


def test_partition_invariants_from_csv(surface_case):
    """The written output mirrors the constructor partition invariants."""
    station = _read_station(surface_case)
    nodes = _read_nodes(surface_case)

    # One element -> one patch row per time step, every patch non-empty
    assert len(station) == N_TIMES
    assert np.all(station.area > 0)

    for time in station.time:
        patch_area = station[station.time == time].area.sum()
        node_area = nodes[nodes.time == time].area.sum()
        assert_allclose(patch_area, PLATE_AREA, rtol=AREA_RTOL)
        assert_allclose(patch_area, node_area, rtol=1e-12)
        assert nodes[nodes.time == time].node.nunique() == N_NODES

    # Every node is assigned to the single patch: its station lies within the
    # patch span and the element station lies inside it (non-empty patch)
    assert nodes.station.min() < ELEMENT_STATION < nodes.station.max()


def test_instrumentation_line(surface_case, default_case):
    """Per-`addSup` instrumentation reports the bounded candidate query (D6).

    The instrumented run emits one pinned line per `distribute()` and appends
    a matching row to `<owner>.surface_distribution.csv`. The candidate count
    must be far below the naive full local-cell scan (`N_NODES *
    N_LOCAL_CELLS`) and the instrumentation must not alter the distribution:
    the patch force is still the analytic element force and the field integral
    still equals the patch total (partition of unity). The no-surface default
    run is never instrumented.
    """
    log = _read_log(surface_case)
    matches = DISTRIBUTION_RE.findall(log)
    assert matches, "no per-addSup instrumentation line in the surface run"

    # The candidate lists are built once, so every call reports the same counts
    counts = set()
    for nodes, candidates, mean, max_candidates, seconds in matches:
        assert int(nodes) == N_NODES
        assert int(candidates) > 0
        assert int(max_candidates) >= 1
        assert float(seconds) >= 0.0
        assert_allclose(
            float(mean), int(candidates) / N_NODES, rtol=1e-4, atol=1e-9
        )
        counts.add((int(candidates), int(max_candidates)))
    assert len(counts) == 1

    candidates, max_candidates = counts.pop()

    # Bounded query: the support stencil (<= 5^3 cells per node here) is far
    # below the naive N_NODES * N_LOCAL_CELLS scan
    naive = N_NODES * N_LOCAL_CELLS
    assert candidates < naive / 10
    assert candidates <= N_NODES * 5**3
    assert max_candidates <= 5**3
    assert max_candidates >= candidates // N_NODES

    # The CSV mirrors the line, with the documented schema
    dist = _read_distribution(surface_case)
    assert list(dist.columns) == DISTRIBUTION_COLUMNS
    assert len(dist) == len(matches)
    assert_array_equal(dist.nodes, N_NODES)
    assert_array_equal(dist.candidates, candidates)
    assert_array_equal(dist.max_candidates, max_candidates)
    assert np.all(np.isfinite(dist.to_numpy(dtype=float)))
    assert np.all(dist.seconds >= 0.0)

    # Instrumentation does not change the distribution result: the surface
    # patch force is still the analytic element force and the field integral
    # still equals the patch total (partition of unity)
    station = _read_station(surface_case)
    assert_allclose(
        _row_vector(station.iloc[0], ["force_x", "force_y", "force_z"]),
        [F_ELEMENT, 0.0, 0.0],
        rtol=1e-12,
        atol=1e-15,
    )
    for time, vec in _read_force_integral(surface_case):
        patch = _row_vector(
            station[station.time == time].iloc[-1],
            ["force_x", "force_y", "force_z"],
        )
        assert_allclose(vec, -patch, rtol=1e-5, atol=1e-12)

    # The no-surface default run is never instrumented
    assert not DISTRIBUTION_RE.search(_read_log(default_case))
    assert not os.path.isfile(os.path.join(default_case, DISTRIBUTION_CSV))


# Rotating-fixture tests (design 8.1): the minimal AFTAL fixture drives the
# rotation lockstep (test_rotation_lockstep) and the D8 surface-moment
# replacement (test_surface_moment_in_torque). The standalone fixture drives
# the MPI total, the Gaussian ablation and the STL failure paths below, and
# system/fvOptions.surface.gaussian backs test_kernel_gaussian_conserves.


def test_rotation_lockstep(aftal_case):
    """The surface nodes follow the rotor rotation exactly (D5).

    ``turbine.csv`` reports the azimuth at both time steps; rotating the first
    step's global node positions by the reported delta must reproduce the
    second step's, and the element positions must rotate by the same delta, so
    the surface is in lockstep with the elements rather than merely close to
    them.
    """
    nodes = _read_aftal_nodes(aftal_case)
    turbine = _read_turbine(aftal_case)
    element = _read_aftal_element(aftal_case)

    # Stable node set and count at both steps
    for time in turbine.time:
        assert (nodes.time == time).sum() == AFTAL_N_NODES
    assert sorted(nodes.node.unique()) == list(range(AFTAL_N_NODES))
    assert len(element) == len(turbine)
    assert np.all(np.isfinite(element[["x", "y", "z"]].to_numpy(dtype=float)))

    t1, t2 = turbine.time.iloc[0], turbine.time.iloc[-1]
    delta = np.deg2rad(
        float(turbine[turbine.time == t2].iloc[0].angle_deg)
        - float(turbine[turbine.time == t1].iloc[0].angle_deg)
    )

    # The reported azimuth delta is the rotor advance over the step
    assert_allclose(delta, OMEGA * (t2 - t1), rtol=1e-4)

    basis = _body_basis([0.0, 0.0, 1.0])
    X1 = _to_global(
        nodes[nodes.time == t1].sort_values("node")[["x", "y", "z"]], basis
    )
    X2 = _to_global(
        nodes[nodes.time == t2].sort_values("node")[["x", "y", "z"]], basis
    )
    assert not np.allclose(X1, X2)

    # Tolerance: turbine.csv writes the azimuth and the element CSV the
    # positions with six significant digits, so the rounded delta dominates
    # the error (about 2e-7 m here), far below this bound
    rotation = _axis_rotation(delta)
    assert_allclose(X2, (rotation @ X1.T).T, rtol=0, atol=1e-5)

    # The elements rotate by the same delta: the surface follows them exactly
    Xe1 = element[element.time == t1].iloc[0][["x", "y", "z"]].to_numpy(float)
    Xe2 = element[element.time == t2].iloc[0][["x", "y", "z"]].to_numpy(float)
    assert_allclose(Xe2, rotation @ Xe1, rtol=0, atol=1e-5)
    assert_allclose(
        element[element.time == t1].iloc[0].root_dist, 0.5, rtol=0, atol=1e-12
    )


def test_surface_moment_in_torque(aftal_case):
    """The reported torque is the surface node moment (D8), not the element one.

    ``turbine.csv`` writes ``ct``; the fixture torque is
    ``ct*0.5*frontalArea*rotorRadius*|Vinf|^2``. The per-node CSV forces are
    per unit density (rho 1.0), the same quantities ``moment()`` uses, so the
    rotor-axis projection of ``sum X_i x F_i`` must reconstruct it at
    ``rtol 1e-6``. The element-position moment (what the delivered ``moment()``
    would report) is more than 20% away: the D8 replacement convention is
    asserted, not merely that a moment exists.
    """
    nodes = _read_aftal_nodes(aftal_case)
    turbine = _read_turbine(aftal_case)
    element = _read_aftal_element(aftal_case)

    basis = _body_basis([0.0, 0.0, 1.0])

    for time in turbine.time:
        d = nodes[nodes.time == time].sort_values("node")
        X = _to_global(d[["x", "y", "z"]], basis)
        F = _to_global(d[["fx", "fy", "fz"]], basis)
        moment = np.sum(np.cross(X, F), axis=0)

        torque = float(turbine[turbine.time == time].iloc[0].ct) * (
            TORQUE_DENOMINATOR
        )
        assert_allclose(moment @ ROTOR_AXIS, torque, rtol=1e-6, atol=0.0)

        # The moment of the same force applied at the element position is far
        # outside the band: the surface moment is the load actually applied
        e = element[element.time == time].iloc[0]
        element_moment = np.cross(
            [float(e.x), float(e.y), float(e.z)], F.sum(axis=0)
        )
        assert not np.isclose(element_moment @ ROTOR_AXIS, torque, rtol=1e-6)
        assert abs(element_moment @ ROTOR_AXIS - torque) > 0.05 * abs(torque)


def test_parallel_total_preserved(surface_case, tmp_path):
    """Two ranks: the distributed total is preserved across the split.

    The standalone fixture's processor boundary is x = 0.5 (``n (2 1 1)``),
    i.e. the plate plane, so both nodes' 5^3-cell supports straddle it: each
    rank sees only its local candidates and the global total must equal the
    serial run's.
    """
    case_dir = str(tmp_path / "parallel")
    _copy_case(CASE_DIR, case_dir)
    _run_case(case_dir, "-surface", "-parallel")

    log = _read_log(case_dir)
    assert "Finalising parallel run" in log
    assert os.path.isdir(os.path.join(case_dir, "processor0"))
    assert os.path.isdir(os.path.join(case_dir, "processor1"))

    decomp = open(os.path.join(case_dir, "system", "decomposeParDict")).read()
    assert "numberOfSubdomains 2" in decomp
    assert "n               (2 1 1)" in decomp

    serial_station = _read_station(surface_case)
    parallel_station = _read_station(case_dir)
    assert list(parallel_station.columns) == STATION_COLUMNS

    for time in serial_station.time:
        serial_patch = _row_vector(
            serial_station[serial_station.time == time].iloc[-1],
            ["force_x", "force_y", "force_z"],
        )
        parallel_patch = _row_vector(
            parallel_station[parallel_station.time == time].iloc[-1],
            ["force_x", "force_y", "force_z"],
        )
        assert_allclose(parallel_patch, serial_patch, rtol=1e-6, atol=1e-9)

    # The field integral equals the serial total and the parallel patch total
    serial_integral = _read_force_integral(surface_case)
    parallel_integral = _read_force_integral(case_dir)
    assert len(parallel_integral) == len(serial_integral)
    for (t_serial, v_serial), (t_parallel, v_parallel) in zip(
        serial_integral, parallel_integral
    ):
        assert t_parallel == t_serial
        assert_allclose(v_parallel[0], v_serial[0], rtol=1e-6, atol=1e-12)
        patch = _row_vector(
            parallel_station[parallel_station.time == t_parallel].iloc[-1],
            ["force_x", "force_y", "force_z"],
        )
        assert_allclose(v_parallel, -patch, rtol=1e-5, atol=1e-9)


def test_kernel_gaussian_conserves(gaussian_case, surface_case):
    """The Gaussian ablation conserves the total (width change only, D7).

    The truncated-Gaussian width is the no-mesh ASM epsilon; its discrete
    volume sum is the pinned ``GAUSSIAN_KERNEL_SUM``, so the field integral is
    that fraction of the patch total while the per-station force is still the
    analytic element force. The cosine run keeps its own (exact) sum: the
    ablation changes the width, not the load.
    """
    log = _read_log(gaussian_case)
    assert "kernel gaussian" in log
    assert "meshFactor 1" in log

    station = _read_station(gaussian_case)
    assert_allclose(
        _row_vector(station.iloc[0], ["force_x", "force_y", "force_z"]),
        [F_ELEMENT, 0.0, 0.0],
        rtol=1e-12,
        atol=1e-15,
    )
    assert np.all(np.isfinite(station.to_numpy(dtype=float)))

    integral = _read_force_integral(gaussian_case)
    assert len(integral) == N_TIMES
    for time, vec in integral:
        patch = _row_vector(
            station[station.time == time].iloc[-1],
            ["force_x", "force_y", "force_z"],
        )
        assert_allclose(vec, -patch, rtol=1e-3, atol=1e-12)
        assert_allclose(-vec[0] / patch[0], GAUSSIAN_KERNEL_SUM, rtol=1e-4)

    # The ablation is a different width from the cosine default: the
    # distributed total differs while the patch total does not (the
    # comparison isolates geometry from the kernel/width confound)
    cosine_first = _read_force_integral(surface_case)[0][1]
    gaussian_first = integral[0][1]
    assert not np.isclose(gaussian_first[0], cosine_first[0], rtol=1e-4)
    assert_allclose(
        _read_station(surface_case).iloc[0].force_x,
        station.iloc[0].force_x,
        rtol=1e-12,
    )


def _install_surface_variant(case_dir):
    """Install ``system/fvOptions.surface`` over the default fvOptions."""
    shutil.copyfile(
        os.path.join(case_dir, "system", "fvOptions.surface"),
        os.path.join(case_dir, "system", "fvOptions"),
    )


def _run_with_surface_geometry(case_dir, geometry):
    """Point the surface at `geometry`, run the solver, return (code, log)."""
    fv_options = os.path.join(case_dir, "system", "fvOptions")
    text = open(fv_options).read()
    assert "geometry/plate.stl" in text
    with open(fv_options, "w") as fv:
        fv.write(text.replace("geometry/plate.stl", geometry))

    subprocess.check_call(["cp", "-rf", "0.org", "0"], cwd=case_dir)
    for cmd in ("blockMesh", "topoSet"):
        subprocess.check_call(cmd, cwd=case_dir, stdout=subprocess.DEVNULL)

    log_path = os.path.join(case_dir, "log.pimpleFoam")
    with open(log_path, "w") as log:
        proc = subprocess.run(
            ["pimpleFoam"],
            cwd=case_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    return proc.returncode, open(log_path).read()


def test_missing_empty_corrupt_stl_aborts(tmp_path):
    """Missing, empty and unparseable surfaceGeometry abort (FatalError)."""
    # Missing file
    case_dir = _copy_case(CASE_DIR, str(tmp_path / "missing"))
    _install_surface_variant(case_dir)
    code, log = _run_with_surface_geometry(
        case_dir, "geometry/doesNotExist.stl"
    )
    assert code != 0
    assert "FOAM FATAL ERROR" in log
    assert "not found" in log

    # Empty file: a valid STL carrying no triangles
    case_dir = _copy_case(CASE_DIR, str(tmp_path / "empty"))
    _install_surface_variant(case_dir)
    code, log = _run_with_surface_geometry(case_dir, "geometry/empty.stl")
    assert code != 0
    assert "FOAM FATAL ERROR" in log
    assert "contains no triangles" in log

    # Unparseable file
    case_dir = _copy_case(CASE_DIR, str(tmp_path / "corrupt"))
    _install_surface_variant(case_dir)
    with open(os.path.join(case_dir, "geometry", "corrupt.stl"), "w") as stl:
        stl.write("this is not an STL file\n")
    code, log = _run_with_surface_geometry(case_dir, "geometry/corrupt.stl")
    assert code != 0
    assert "FOAM FATAL ERROR" in log
