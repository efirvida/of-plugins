# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-Python sanity checks for the committed NREL Phase VI data.

Covers the geometry, the measured experimental anchors, the polar files, the
per-directory provenance and the no-report-artefact rule (spec requirements
D1-D7). These tests never require OpenFOAM.
"""

from __future__ import annotations

import csv
import hashlib
import os
import re
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
TURBINESFOAM = TESTS_DIR.parent
PACKAGE = TURBINESFOAM / "validation" / "phaseVI"
DATA = PACKAGE / "data"
sys.path.insert(0, str(PACKAGE / "tools"))

from case_config import (  # noqa: E402
    assert_tsr_matches_experiment,
    load_config,
    read_blade,
)
from element_data import blade_element_rows  # noqa: E402

CONFIG = PACKAGE / "config" / "case.yaml"
GEOMETRY_SHA256 = (
    "e678d6c52df9c8562361d61c65ad373535ea2351448874213b706bc6824e1576"
)
H_PERFORMANCE_ANCHORS = {
    "7.0": ("h0700000", 5.9458, 789.8572, 1132.0692),
    "25.0": ("h2500000", 11.9507, 1580.4082, 4028.6306),
}
NO_ARTEFACT_SUFFIXES = {".pdf", ".xls", ".xlsx", ".xlsm", ".doc", ".docx"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def polar_rows(text: str) -> list[tuple[float, ...]]:
    rows = []
    for line in text.splitlines():
        match = re.match(r"^\(([-\d.eE+ ]+)\)$", line.strip())
        if match:
            rows.append(tuple(float(value) for value in match.group(1).split()))
    return rows


def polar_blocks(text: str) -> dict[str, list[tuple[float, ...]]]:
    blocks: dict[str, list[tuple[float, ...]]] = {}
    current = None
    for line in text.splitlines():
        header = re.match(r"^(\w+Data)$", line.strip())
        if header:
            current = header.group(1)
            blocks[current] = []
            continue
        match = re.match(r"^\(([-\d.eE+ ]+)\)$", line.strip())
        if current is not None and match:
            blocks[current].append(
                tuple(float(value) for value in match.group(1).split())
            )
    return blocks


def test_geometry_file_hash_and_shape():
    path = DATA / "geometry" / "phaseVI_blade.csv"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == GEOMETRY_SHA256
    rows = read_blade(path)
    assert len(rows) == 26
    assert rows[0]["radius"] == pytest.approx(0.5083)
    assert rows[-1]["radius"] == pytest.approx(5.029)
    for previous, current in zip(rows, rows[1:]):
        assert current["radius"] > previous["radius"]
    for row in rows:
        assert row["chord"] > 0.0
        assert 0.0 < row["mount"] <= 0.5
    assert rows[0]["chord"] == pytest.approx(0.2180)
    assert rows[-1]["chord"] == pytest.approx(0.3551)
    assert rows[0]["mount"] == pytest.approx(0.50)
    assert rows[-1]["mount"] == pytest.approx(0.30)


def test_geometry_mounting_convention():
    cfg = load_config(CONFIG)
    pitch = float(cfg["turbine"]["pitch_deg"])
    geometry = read_blade()
    elements = blade_element_rows(cfg)
    assert len(elements) == len(geometry)
    for row, element in zip(geometry, elements):
        assert element[1] == pytest.approx(row["radius"])
        assert element[3] == pytest.approx(row["chord"])
        assert element[4] == pytest.approx(row["mount"])
        # the element writer converts the report twist with -(twist + pitch)
        assert element[5] == pytest.approx(-(row["twist"] + pitch))
    assert any(row["twist"] > 10.0 for row in geometry)


def test_experimental_anchors():
    rows = {
        row["wind_speed_m_s"]: row
        for row in read_csv(DATA / "experiment" / "sequence_H_performance.csv")
    }
    assert sorted(rows, key=float) == ["7.0", "10.0", "13.0", "15.0", "20.0", "25.0"]
    for speed, (name, power, torque, thrust) in H_PERFORMANCE_ANCHORS.items():
        row = rows[speed]
        assert row["file_name"] == name
        assert float(row["rotpow_kw"]) == pytest.approx(power, abs=1e-4)
        assert float(row["lsstqcor_nm"]) == pytest.approx(torque, abs=1e-4)
        assert float(row["eaeroth_n"]) == pytest.approx(thrust, abs=1e-4)
        assert abs(float(row["yaw_deg"])) <= 0.5
        assert int(row["repetition"]) == 0
    corrected = rows["20.0"]
    assert corrected["file_name"] == "h20m0000"
    assert float(corrected["rotpow_kw"]) == pytest.approx(9.3661, abs=1e-4)


def test_sequence_s_repeat():
    rows = read_csv(DATA / "experiment" / "sequence_S_performance.csv")
    assert len(rows) == 1
    row = rows[0]
    assert row["file_name"] == "s0700000"
    assert float(row["wind_speed_m_s"]) == pytest.approx(7.0)
    assert float(row["rotpow_kw"]) == pytest.approx(6.0304, abs=1e-4)
    assert float(row["lsstqcor_nm"]) == pytest.approx(801.265, abs=1e-4)
    assert float(row["eaeroth_n"]) == pytest.approx(1149.5288, abs=1e-4)


def test_spanwise_stations():
    rows = read_csv(DATA / "experiment" / "sequence_H_spanwise.csv")
    stations = sorted({float(row["r_over_R"]) for row in rows})
    assert stations == [0.3, 0.47, 0.63, 0.8, 0.95]
    for row in rows:
        assert float(row["cn"]) > 0.0
        assert float(row["cm"]) < 0.0
    # the 7 m/s sign gate relies on both measured coefficients being positive
    seven = [row for row in rows if float(row["wind_speed_m_s"]) == 7.0]
    assert len(seven) == 5
    assert all(float(row["ct"]) > 0.0 for row in seven)


def test_tsr_matches_measured_rows():
    cfg = load_config(CONFIG)
    assert_tsr_matches_experiment(cfg)


def test_polar_baseline():
    path = DATA / "polars" / "S809_OSU_Re1M_total.dat"
    text = path.read_text(encoding="utf-8")
    assert "Table A-7" in text
    rows = polar_rows(text)
    assert len(rows) == 59
    alpha = [row[0] for row in rows]
    assert alpha == sorted(alpha)
    assert alpha[0] == pytest.approx(-180.0)
    assert alpha[-1] == pytest.approx(180.0)
    assert all(row[2] >= 0.0 for row in rows)
    by_alpha = {row[0]: row for row in rows}
    assert by_alpha[8.2][1] == pytest.approx(0.90)


def test_polar_multire():
    text = (DATA / "polars" / "S809_multiRe.dat").read_text(encoding="utf-8")
    assert "tableType multiRe;" in text
    re_list = re.search(r"ReList\s*\(([^)]*)\)", text)
    assert re_list is not None
    levels = [float(value) for value in re_list.group(1).split()]
    assert levels == [7.5e5, 1.0e6, 1.25e6, 1.5e6]
    blocks = polar_blocks(text)
    assert set(blocks) == {"clData", "cdData", "cmData"}
    for name, rows in blocks.items():
        # one row per shared alpha station and four Re columns per row
        assert len(rows) == 28, name
        assert all(len(row) == 5 for row in rows)
        alpha = [row[0] for row in rows]
        assert alpha == sorted(alpha)
        assert len(set(alpha)) == len(alpha)
    provenance = (DATA / "polars" / "PROVENANCE.md").read_text(encoding="utf-8")
    assert "S809_Re1M_extended" in provenance
    assert "is **not** used" in provenance


def test_provenance_fields():
    geometry = (DATA / "geometry" / "PROVENANCE.md").read_text(encoding="utf-8")
    polars = (DATA / "polars" / "PROVENANCE.md").read_text(encoding="utf-8")
    experiment = (DATA / "experiment" / "PROVENANCE.md").read_text(encoding="utf-8")
    for text in (geometry, polars, experiment):
        assert "2026-09-19" in text
        assert "8284be8c" in text
        assert "52fd258" in text
    assert "10.2172/15000240" in geometry
    assert GEOMETRY_SHA256 in geometry
    assert "Table A-1" in geometry
    assert "8e1cfad2" in experiment
    assert "ldsmean" in experiment
    assert "WTBARO" in experiment


def test_no_report_artefacts_committed():
    offenders = []
    for root, dirs, files in os.walk(PACKAGE):
        dirs[:] = [
            name
            for name in dirs
            if name not in {"runs", "results", "__pycache__"}
        ]
        for name in files:
            if Path(name).suffix.lower() in NO_ARTEFACT_SUFFIXES:
                offenders.append(str(Path(root) / name))
    assert offenders == []
