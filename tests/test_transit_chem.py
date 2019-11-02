#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `transit_chem` package."""


from click.testing import CliRunner

import transit_chem
from transit_chem import cli


def test_content():
    transit_chem.HarmonicOscillator


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    help_result = runner.invoke(cli.main, ["--help"])
    assert help_result.exit_code == 0
