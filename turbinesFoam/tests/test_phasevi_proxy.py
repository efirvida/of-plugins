# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-Python contracts for the rotational-augmentation proxy harness.

These tests never require OpenFOAM and never submit a job: they exercise the
committed variant matrix, the serial-chain plan, the control gate and the
fail-loud criteria of
`validation/phaseVI/scripts/proxyRotationalAugmentation.py` (spec
`phasevi-proxy-verification`, design sections 7-8).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
TURBINESFOAM = TESTS_DIR.parent
PACKAGE = TURBINESFOAM / "validation" / "phaseVI"
sys.path.insert(0, str(PACKAGE / "scripts"))
sys.path.insert(0, str(PACKAGE / "tools"))

import proxyRotationalAugmentation as proxy  # noqa: E402
from case_config import load_config  # noqa: E402

CONFIG = PACKAGE / "config" / "case.yaml"


@pytest.fixture(scope="module")
def cfg():
    return load_config(CONFIG)


def synth_metrics(
    control_cp=-0.0459,
    control_ct=-0.0157,
    control_c_ref_t=-0.063,
    augment_cp=0.01,
    augment_ct=0.01,
    augment_c_ref_t=0.03,
    control_power_kw=5.17,
    ablation_power_kw=5.57,
):
    """A complete synthetic metric set for `check_gates`."""
    return {
        "13": {
            "control": {
                "cp": control_cp,
                "ct": control_ct,
                "midspan_c_ref_t": control_c_ref_t,
            },
            "augmentation-on": {
                "cp": augment_cp,
                "ct": augment_ct,
                "midspan_c_ref_t": augment_c_ref_t,
                "torque_nm": augment_ct * 1.0,
            },
            "augmentation-on-root-off": {},
        },
        "7": {
            "control": {"power_kw": control_power_kw},
            "augmentation-on": {},
            "augmentation-on-root-off": {"power_kw": ablation_power_kw},
        },
    }


def synth_measured():
    return {"7": {"rotpow_kw": 5.9458}}


def test_variant_matrix(cfg):
    """Exactly three variants, differing only in the two toggles."""
    assert [variant["name"] for variant in proxy.VARIANTS] == [
        "control",
        "augmentation-on",
        "augmentation-on-root-off",
    ]
    proxy.assert_variants_share_keys(cfg)

    matrix = proxy.variant_matrix(cfg)
    assert len(matrix) == len(proxy.SPEEDS) * 3
    for speed in proxy.SPEEDS:
        stripped = {
            variant["name"]: proxy.toggle_stripped(matrix[(speed, variant["name"])])
            for variant in proxy.VARIANTS
        }
        assert stripped["control"] == stripped["augmentation-on"]
        assert stripped["control"] == stripped["augmentation-on-root-off"]
        assert "active off;" in matrix[(speed, "control")]
        assert "active on;" in matrix[(speed, "augmentation-on")]
        assert "rootEffects on;" in matrix[(speed, "augmentation-on")]
        assert "rootEffects off;" in matrix[(speed, "augmentation-on-root-off")]


def test_check_mode_reports_proxy_and_submits_nothing(cfg, monkeypatch, capsys):
    """`--check` prints the 7/13 m/s, D/32, 0.25-rev configuration."""
    calls = []

    def guard(*args, **kwargs):
        calls.append((args, kwargs))
        raise AssertionError(f"check mode ran a subprocess: {args}")

    monkeypatch.setattr(proxy.subprocess, "run", guard)
    assert proxy.main(["--check"]) == 0
    captured = capsys.readouterr()
    payload = captured.out.split("check:", 1)[0]
    config = json.loads(payload)
    assert config["speeds_m_s"] == ["7", "13"]
    assert config["mesh_cell_size_D"] == 32
    assert config["revolutions"] == 0.25
    assert config["queue"] == "sequana_cpu_dev"
    assert config["production_submission"] is False
    assert [variant["name"] for variant in config["variants"]] == [
        "control",
        "augmentation-on",
        "augmentation-on-root-off",
    ]
    assert calls == []


def test_serial_chain():
    """The variants are chained serially, never submitted as an array."""
    plan = proxy.execution_plan()
    assert plan["mode"] == "serial"
    assert plan["array"] is False
    assert plan["max_submit"] == 1
    assert plan["steps"] == [
        {"speed": speed, "variant": variant["name"]}
        for speed in proxy.SPEEDS
        for variant in proxy.VARIANTS
    ]

    command = proxy.build_submit_command(
        proxy.SCRIPTS / "proxyRotationalAugmentation.sh", proxy.DEV_QUEUE
    )
    assert command[0] == "sbatch"
    assert f"--partition={proxy.DEV_QUEUE}" in command
    assert not any(part == "--array" for part in command)
    assert not any("--array" in part for part in command)


def test_fail_loud(cfg):
    """A control `cp` outside tolerance fails and names the control gate."""
    failures, report = proxy.check_gates(
        synth_metrics(control_cp=-0.06), synth_measured()
    )
    assert failures
    assert any("control gate" in failure for failure in failures)
    assert report["control_gate"]["pass"] is False

    # A negative augmentation cp fails the primary signal and names the variant.
    failures, report = proxy.check_gates(
        synth_metrics(augment_cp=-0.02), synth_measured()
    )
    assert any("U13 augmentation-on primary" in failure for failure in failures)
    assert report["primary"]["U13_augmentation_on"]["cp_positive"] is False

    # A non-positive mid-span c_ref_t fails the secondary signal.
    failures, _ = proxy.check_gates(
        synth_metrics(augment_c_ref_t=-0.01), synth_measured()
    )
    assert any("secondary signal" in failure for failure in failures)

    # A non-shrinking U7 deficit fails the primary signal.
    failures, _ = proxy.check_gates(
        synth_metrics(control_power_kw=5.57, ablation_power_kw=5.17),
        synth_measured(),
    )
    assert any("U7 augmentation-on + root-off" in failure for failure in failures)

    # A fully met synthetic set passes with no failures.
    failures, report = proxy.check_gates(synth_metrics(), synth_measured())
    assert failures == []
    assert report["pass"] is True
    assert report["primary"]["U13_augmentation_on"]["cp_positive"] is True
    assert report["secondary"]["U13_augmentation_on"]["pass"] is True
    assert report["primary"]["U7_root_off_ablation"]["deficit_shrunk"] is True


def test_no_production_submission(cfg):
    """No code path submits to a production queue; no array is ever built."""
    for queue in ("sequana_cpu_long", "production", "prod-queue", "long"):
        with pytest.raises(ValueError):
            proxy.build_submit_command(Path("job.sh"), queue)
    # The one accepted queue is the authorized development queue.
    assert proxy.build_submit_command(Path("job.sh"), proxy.DEV_QUEUE)[0] == "sbatch"

    source = (proxy.SCRIPTS / "proxyRotationalAugmentation.py").read_text()
    assert "--array" not in source
    assert "sequana_cpu_long" not in source
    assert "production.slurm" not in source


def test_prepare_isolates_variant_toggles(tmp_path, cfg, monkeypatch):
    """Each prepared copy keeps its own fvOptions despite the hardlink base."""
    def fake_base(work_root, speed, ranks, force):
        base = proxy.base_run_dir(work_root, speed)
        (base / "system").mkdir(parents=True)
        # The base twin is hardlinked into every variant; a shared inode would
        # make the last write win.
        (base / "system" / "fvOptions").write_text("base\n")
        (base / "0").mkdir()
        (base / "constant" / "polyMesh").mkdir(parents=True)
        (base / "constant" / "polyMesh" / "points").write_text("mesh\n")
        return base

    monkeypatch.setattr(proxy, "prepare_base", fake_base)
    prepared = proxy.prepare(tmp_path, cfg, proxy.RUN_RANKS)
    assert len(prepared) == len(proxy.SPEEDS) * len(proxy.VARIANTS)
    for speed in proxy.SPEEDS:
        for variant in proxy.VARIANTS:
            vdir = proxy.variant_dir(tmp_path, speed, variant)
            text = (vdir / "system" / "fvOptions").read_text()
            proxy.assert_rendered_toggle(text, variant)
            # The base twin must stay untouched (the inode was not shared).
            base = proxy.base_run_dir(tmp_path, speed)
            assert (base / "system" / "fvOptions").read_text() == "base\n"


def test_measured_anchors():
    """The measured U7/U13 anchors used by the gates are read from the case."""
    measured = proxy.load_measured()
    assert measured["7"]["rotpow_kw"] == pytest.approx(5.9458)
    assert measured["13"]["rotpow_kw"] == pytest.approx(9.7936)
    assert measured["7"]["rho_kg_m3"] > 0.0
