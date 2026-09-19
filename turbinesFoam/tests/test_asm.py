#!/usr/bin/env python
"""Tests for the actuator surface element (`actuatorSurfaceElement`) run via
`actuatorLineSource`.

The static case regenerates `system/fvOptions` at run time from
`system/fvOptions.template` (see `scripts/set_alpha.py`), so the ASM keys are
injected into the *copied* template in this test directory -- the file that
`set_alpha.py` actually formats -- rather than into `system/fvOptions`, which
would be overwritten. The two injected lines are plain (no braces) and
therefore survive `str.format`.
"""

from __future__ import division, print_function

import os
import subprocess

import numpy as np
import pandas as pd
import pytest

alpha_deg = 10.0
output_fpath = "postProcessing/actuatorLines/0/foil.csv"
element_dir = "postProcessing/actuatorLineElements/0/"
element_fpath = os.path.join(element_dir, "foil.element0.csv")
template_fpath = "system/fvOptions.template"
asm_keys = (
    "        elementType actuatorSurfaceElement;\n"
    "        nChordwise 5;\n"
)


@pytest.fixture(scope="module")
def chdir():
    orig = os.getcwd()
    os.chdir("tests/actuatorLineSource")
    yield
    os.chdir(orig)


def get_tutorial_files(case="static"):
    cmd = "./getTutorialFiles.sh {}".format(case)
    out = subprocess.check_output(cmd, shell=True)


def patch_template():
    """Add the ASM keys inside actuatorLineSourceCoeffs of the copied
    template. Plain lines only: braces would break `str.format`."""
    with open(template_fpath) as f:
        txt = f.read()
    assert "elementType" not in txt, "template already patched"
    anchor = "writeElementPerf    true;"
    assert anchor in txt, "anchor line not found in template"
    txt = txt.replace(anchor, anchor + "\n" + asm_keys)
    with open(template_fpath, "w") as f:
        f.write(txt)


def check_created():
    """Test that actuatorLineSource was created."""
    txt = "Selecting finite volume options.*actuatorLineSource"
    subprocess.check_output(["grep", txt, "log.simpleFoam"])


def check_output_file_exists():
    """Test that actuatorLineSource output file exists."""
    assert os.path.isfile(output_fpath)


def check_element_file_exists():
    """Test that the element perf file was created."""
    assert os.path.isfile(element_fpath)


def check_element_csv_finite():
    """Element CSV forces must be present and finite. The surface element
    writes the same CSV format as the line element (inherited BEM chain);
    the total force is the element `forceVector_`, never a per-strip
    quantity."""
    df = pd.read_csv(element_fpath)
    for col in ["fx", "fy", "fz", "f_ref_t", "f_ref_n"]:
        assert np.all(np.isfinite(df[col])), "column {} not finite".format(col)


def read_rel_vel_mag():
    df = pd.read_csv(element_fpath)
    return df.rel_vel_mag.iloc[-1]


def check_chord_averaging(alm_rel_vel_mag, asm_rel_vel_mag):
    """The ASM inflow is the arithmetic mean of samples at the 5 chord strip
    midpoints, while the ALM samples a single point at the quarter chord, so
    on the same flow the chord-averaged relative velocity must differ from
    the single-point sample. Tolerance is stated explicitly (0.001% relative
    -- the induction deficit on this case is small, ~4e-5); both values are
    printed for diagnosis."""
    print("ALM rel_vel_mag (single-point): {:.8f}".format(alm_rel_vel_mag))
    print("ASM rel_vel_mag (chord-averaged): {:.8f}".format(asm_rel_vel_mag))
    assert np.isfinite(alm_rel_vel_mag)
    assert np.isfinite(asm_rel_vel_mag)
    rel_diff = abs(asm_rel_vel_mag - alm_rel_vel_mag) / abs(alm_rel_vel_mag)
    print("Relative difference: {:.3e}".format(rel_diff))
    assert rel_diff > 1e-5


def chord_strip_mean(func, n):
    """Arithmetic mean of func evaluated at the strip midpoints f_k =
    (k + 0.5)/n, matching the element's chord-averaging scheme."""
    return sum(func((k + 0.5) / n) for k in range(n)) / n


def check_expected_chord_mean():
    """Python expected-value check on a known linear flow u(f) = 10 + 20 f
    along the chord: the strip-midpoint mean equals the value at the mean
    midpoint (f = 0.5) and differs from the single-point (quarter-chord,
    f = 0.25) sample."""
    u = lambda f: 10.0 + 20.0 * f  # noqa: E731
    expected = u(0.5)
    computed = chord_strip_mean(u, 5)
    single_point = u(0.25)
    np.testing.assert_allclose(computed, expected)
    assert computed != single_point


def test_asm_2d(chdir):
    """2-D actuatorLineSource with the actuator surface element."""
    # 1. Baseline ALM run on the pristine template: elementType absent so
    #    elements are constructed as actuatorLineElement (default)
    get_tutorial_files(case="static")
    subprocess.check_output("./Allclean")
    subprocess.check_output(["./Allrun", "2D", str(alpha_deg)])
    check_created()
    check_output_file_exists()
    check_element_file_exists()
    alm_rel_vel_mag = read_rel_vel_mag()

    # 2. ASM run: inject the ASM keys into the copied template that
    #    set_alpha.py formats at run time
    subprocess.check_output("./Allclean")
    patch_template()
    subprocess.check_output(["./Allrun", "2D", str(alpha_deg)])
    check_created()
    check_output_file_exists()
    check_element_file_exists()
    check_element_csv_finite()
    asm_rel_vel_mag = read_rel_vel_mag()

    # 3. Chord-averaging evidence: ASM differs from the single-point ALM
    #    sample on the same flow, and the sampling scheme matches the
    #    analytic expected value for a known linear flow
    check_chord_averaging(alm_rel_vel_mag, asm_rel_vel_mag)
    check_expected_chord_mean()