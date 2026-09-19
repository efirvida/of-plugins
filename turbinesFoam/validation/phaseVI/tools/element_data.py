#!/usr/bin/env python3
"""Blade and hub element-data generation for the NREL Phase VI case.

turbinesFoam element data rows are
`(axialDistance radius azimuth chord chordMount pitch)` for blades and
`(axialDistance height diameter)` for the hub. The mounting angle uses the
NREL-to-turbinesFoam convention `pitch = -(twist + pitch_deg)`; the NREL
report defines positive twist/pitch towards feather while the axial-flow
element has the opposite sign.

Adapted from `mttbrbr/single-actuator-line` (main 8284be8c..., GPL-3.0-or-later).
"""

from __future__ import annotations

from typing import Any

from case_config import (
    ROOT_CUTOUT_RADIUS,
    ROOT_CYLINDER_END,
    read_blade,
    render_cell_list,
)

# Hub: a lumped cylinder spanning +-root cutout along the rotor axis; the
# diameter is twice the blade root station (TP-500-29955 has no hub table).
HUB_DIAMETER = 2.0 * ROOT_CUTOUT_RADIUS


def blade_element_rows(cfg: dict[str, Any], blade: list[dict[str, float]] | None = None) -> list[list[float]]:
    """Element data rows for one blade, in turbinesFoam column order."""
    pitch = float(cfg["turbine"]["pitch_deg"])
    stations = blade if blade is not None else read_blade()
    rows = []
    for station in stations:
        # NREL twist/pitch is positive towards feather; the axial-flow element
        # convention is the opposite sign.
        turbine_pitch = -(station["twist"] + pitch)
        rows.append(
            [
                0.0,
                station["radius"],
                0.0,
                station["chord"],
                station["mount"],
                turbine_pitch,
            ]
        )
    return rows


def blade_element_profiles(cfg: dict[str, Any], blade: list[dict[str, float]] | None = None) -> list[str]:
    """`cylinder` inboard of the root cutout, `S809` outboard.

    The count of cylinder elements follows the element interpolation used by
    `axialFlowTurbineALSource`, matching the shared scaffold.
    """
    stations = blade if blade is not None else read_blade()
    n_elements = int(cfg["turbine"]["n_elements"])
    radius = float(cfg["turbine"]["radius"])
    span = radius - ROOT_CUTOUT_RADIUS
    n_root = max(1, round(n_elements * (ROOT_CYLINDER_END - ROOT_CUTOUT_RADIUS) / span))
    n_root = min(n_root, n_elements)
    del stations
    return ["cylinder"] * n_root + ["S809"] * (n_elements - n_root)


def hub_element_rows() -> list[list[float]]:
    return [
        [0.0, ROOT_CUTOUT_RADIUS, HUB_DIAMETER],
        [0.0, -ROOT_CUTOUT_RADIUS, HUB_DIAMETER],
    ]


def render_blade_element_data(cfg: dict[str, Any], indent: str = "                    ") -> str:
    return render_cell_list(blade_element_rows(cfg), indent=indent)


def render_hub_element_data(indent: str = "                ") -> str:
    return render_cell_list(hub_element_rows(), indent=indent)


def element_data_summary(cfg: dict[str, Any]) -> dict[str, float | int | str]:
    """Values printed by `scripts/makeElementData.py` and asserted by tests."""
    rows = blade_element_rows(cfg)
    return {
        "stations": len(rows),
        "root_radius": rows[0][1],
        "tip_radius": rows[-1][1],
        "hub_diameter": HUB_DIAMETER,
        "root_mounting_deg": rows[0][5],
        "tip_mounting_deg": rows[-1][5],
    }
