#!/usr/bin/env python
"""Tests for `axialFlowTurbineALSource` with the actuator surface element
(`tutorials/axialFlowTurbineASM`)."""

from __future__ import division, print_function

import os
import subprocess

import pytest

element_dir = "postProcessing/actuatorLineElements/0/"
turbine_dir = "postProcessing/turbines/0/"


@pytest.fixture(scope="module")
def setup():
    orig = os.getcwd()
    os.chdir("tests/axialFlowTurbineASMSource")
    subprocess.run("./getTutorialFiles.sh", shell=True)
    yield
    os.chdir(orig)


def check_created():
    """Test that axialFlowTurbineALSource was created."""
    txt = "Selecting finite volume options.*axialFlowTurbineALSource"
    subprocess.check_output(["grep", txt, "log.pimpleFoam"])


def check_turbine_file_exists():
    """Test that the turbine perf file was created."""
    assert os.path.isfile(os.path.join(turbine_dir, "turbine.csv"))


def check_element_file_exists():
    """Test that the element perf file was created."""
    assert os.path.isfile(
        os.path.join(element_dir, "turbine.blade1.element0.csv")
    )


def check_surface_element_ran():
    """Test that the actuator surface element actually ran: its debug output
    (enabled via debugSwitches) reports the chord-averaged inflow."""
    subprocess.check_output(
        ["grep", "chord-averaged inflow velocity", "log.pimpleFoam"]
    )


def all_checks():
    """All checks to run after simulation."""
    check_created()
    check_turbine_file_exists()
    check_element_file_exists()
    check_surface_element_ran()


def test_serial(setup):
    """Test axialFlowTurbineALSource with the surface element in serial."""
    output_clean = subprocess.check_output("./Allclean")
    try:
        output_run = subprocess.check_output("./Allrun")
    except subprocess.CalledProcessError:
        print(
            subprocess.check_output(
                ["tail", "-n", "200", "log.pimpleFoam"]
            ).decode()
        )
        assert False
    all_checks()
    log_end = subprocess.check_output(["tail", "log.pimpleFoam"]).decode()
    print(log_end)
    assert log_end.split()[-1] == "End"


def test_parallel(setup):
    """Test axialFlowTurbineALSource with the surface element in parallel."""
    output_clean = subprocess.check_output("./Allclean")
    try:
        output_run = subprocess.check_output(["./Allrun", "-parallel"])
    except subprocess.CalledProcessError:
        print(
            subprocess.check_output(
                ["tail", "-n", "200", "log.pimpleFoam"]
            ).decode()
        )
        assert False
    all_checks()
    log_end = subprocess.check_output(["tail", "log.pimpleFoam"]).decode()
    print(log_end)
    assert "Finalising parallel run" in log_end