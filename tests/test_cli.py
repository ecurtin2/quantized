from click.testing import CliRunner

from quantized.cli import cli

import pytest


@pytest.mark.slow
def test_integration_test_via_cli(tmp_path):
    runner = CliRunner()

    runner.invoke(
        cli, ["run-oned", "--cutoff-energy=1", "--until-pnot=0.98", f"--output-path={tmp_path}"]
    )
