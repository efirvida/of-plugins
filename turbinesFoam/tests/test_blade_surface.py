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
rotation/moment tests (W1c2b): one blade, two geometry rows (one element), no
hub/tower, ``writeNodePerf true`` and ``rho 1.0``, with a multi-triangle plate
at the element station so the nodes are offset from the element position.

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
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

CASE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "bladeSurface"
)

STATION_CSV = os.path.join("postProcessing", "bladeSurface", "blade.surface.csv")
NODE_CSV = os.path.join(
    "postProcessing", "bladeSurface", "blade.surface_nodes.csv"
)
ELEMENT_CSV = os.path.join(
    "postProcessing", "actuatorLineElements", "0", "blade.element0.csv"
)

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
ELEMENT_STATION = 0.5
ELEMENT_ROOT_DIST = 0.5

# Triangle centroids in the canonical generation frame (face order of the STL)
NODE_STATION = {0: 0.5333333333333333, 1: 0.4666666666666667}
NODE_CHORD_FRACTION = {0: 0.6666666666666666, 1: 0.3333333333333333}

# Discrete sum of the delivered element strip projection in this fixture
# (truncated Gaussian, unclipped support): the default-path integral is this
# fraction of the element total.
STRIP_KERNEL_SUM = 0.999783


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


def _read_element(case_dir):
    return pd.read_csv(os.path.join(case_dir, ELEMENT_CSV))


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
    assert_allclose(station.area, PLATE_AREA, rtol=1e-12)
    assert np.all(station.area > 0)

    nodes = _read_nodes(surface_case)
    assert list(nodes.columns) == NODE_COLUMNS
    assert len(nodes) == N_TIMES * N_NODES
    assert np.all(np.isfinite(nodes.to_numpy(dtype=float)))

    for time in station.time:
        per_time = nodes[nodes.time == time]
        assert sorted(per_time.node) == [0, 1]
        assert_allclose(per_time.area, NODE_AREA, rtol=1e-12)

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
        assert_allclose(patch_area, PLATE_AREA, rtol=1e-12)
        assert_allclose(patch_area, node_area, rtol=1e-12)
        assert nodes[nodes.time == time].node.nunique() == N_NODES

    # Every node is assigned to the single patch: its station lies within the
    # patch span and the element station lies inside it (non-empty patch)
    assert nodes.station.min() < ELEMENT_STATION < nodes.station.max()


# W1c2b owns the remaining W1.7 tests (design 8.1):
#   - test_rotation_lockstep (tests/bladeSurfaceAFTAL, azimuth delta in
#     turbine.csv vs the per-node positions at two time steps)
#   - test_surface_moment_in_torque (rtol 1e-6, element-position moment must
#     fall outside the band: the D8 replacement, not merely "a moment exists")
#   - test_parallel_total_preserved (mpirun -np 2, node support straddling the
#     x = 0.5 processor boundary)
#   - test_kernel_gaussian_conserves (system/fvOptions.surface.gaussian)
#   - test_missing_empty_corrupt_stl_aborts (geometry/empty.stl + generated
#     corrupt file, mirroring tests/test_nacelle.py)
