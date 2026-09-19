#!/usr/bin/env python
"""Compare the actuator line model (ALM) and actuator surface model (ASM)
tutorial results side by side.

Reads the turbine-level CSV (``postProcessing/turbines/0/turbine.csv``) and
the element-level CSVs (``postProcessing/actuatorLineElements/0/*.csv``) from
two run directories and prints mean Cp/CT/TSR plus spanwise tables for both
models.

Usage::

    python compareALMvsASM.py --alm ../axialFlowTurbineAL --asm .

If either run directory is missing or incomplete, the missing input is
reported and the script exits nonzero instead of fabricating a comparison.
"""

from __future__ import division, print_function

import argparse
import os
import sys

import numpy as np
import pandas as pd

TURBINE_CSV = os.path.join("postProcessing", "turbines", "0", "turbine.csv")
ELEMENT_DIR = os.path.join("postProcessing", "actuatorLineElements", "0")
SPANWISE_COLS = [
    "root_dist", "cl", "alpha_deg", "rel_vel_mag",
    "c_ref_t", "c_ref_n", "f_ref_t", "f_ref_n",
]


def load_run(run_dir, label):
    """Load turbine and element data from a run directory, or report and
    exit nonzero when the run is missing or incomplete."""
    turbine_path = os.path.join(run_dir, TURBINE_CSV)
    element_dir = os.path.join(run_dir, ELEMENT_DIR)
    missing = []
    if not os.path.isfile(turbine_path):
        missing.append(turbine_path)
    if not os.path.isdir(element_dir) or not os.listdir(element_dir):
        missing.append(element_dir + "/*.csv")
    if missing:
        print("Missing input for {} run:".format(label), file=sys.stderr)
        for path in missing:
            print("  {}".format(path), file=sys.stderr)
        print("Run the tutorial case first (./Allrun) and retry.",
              file=sys.stderr)
        sys.exit(1)

    turbine = pd.read_csv(turbine_path).drop_duplicates("time", keep="last")
    elements = []
    for fname in sorted(os.listdir(element_dir)):
        if not fname.endswith(".csv"):
            continue
        df = pd.read_csv(os.path.join(element_dir, fname))
        elements.append(df.iloc[-1])
    spanwise = pd.DataFrame(elements)
    return turbine, spanwise


def print_summary(turbine_alm, turbine_asm):
    """Print mean Cp/CT/TSR side by side."""
    rows = []
    for label, turbine in (("ALM", turbine_alm), ("ASM", turbine_asm)):
        rows.append(
            (
                label,
                turbine.tsr.mean(),
                turbine.cp.mean(),
                turbine.ct.mean(),
            )
        )
    header = "{:<4} {:>10} {:>10} {:>10}".format(
        "Model", "mean TSR", "mean Cp", "mean CT"
    )
    print(header)
    print("-" * len(header))
    for label, tsr, cp, ct in rows:
        print("{:<4} {:>10.4f} {:>10.4f} {:>10.4f}".format(
            label, tsr, cp, ct
        ))


def print_spanwise(span_alm, span_asm):
    """Print spanwise tables side by side (last written time step)."""
    print("")
    print("Spanwise distributions (last time step):")
    for col in SPANWISE_COLS:
        if col not in span_alm.columns or col not in span_asm.columns:
            continue
        header = "{:<9} {:>12} {:>12}".format("r/R", "ALM " + col,
                                              "ASM " + col)
        print(header)
        print("-" * len(header))
        n = min(len(span_alm), len(span_asm))
        for i in range(n):
            print("{:<9.4f} {:>12.6f} {:>12.6f}".format(
                span_alm["root_dist"].iloc[i],
                span_alm[col].iloc[i],
                span_asm[col].iloc[i]
            ))
        print("")


def main():
    parser = argparse.ArgumentParser(
        description="Compare ALM and ASM tutorial results side by side."
    )
    parser.add_argument(
        "--alm", default="../axialFlowTurbineAL",
        help="actuator line model tutorial run directory "
             "(default: ../axialFlowTurbineAL)",
    )
    parser.add_argument(
        "--asm", default=".",
        help="actuator surface model tutorial run directory (default: .)",
    )
    args = parser.parse_args()

    turbine_alm, span_alm = load_run(args.alm, "ALM")
    turbine_asm, span_asm = load_run(args.asm, "ASM")

    print_summary(turbine_alm, turbine_asm)
    print_spanwise(span_alm, span_asm)


if __name__ == "__main__":
    main()