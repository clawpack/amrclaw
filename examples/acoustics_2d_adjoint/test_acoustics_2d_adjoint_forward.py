#!/usr/bin/env python
"""
Regression tests for 2D acoustics with adjoint flagging.
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test
from clawpack.clawutil.test import run_example_for_test

def configure_2d_adjoint(runner):
    clawdata = runner.rundata.clawdata
    gaugedata = runner.rundata.gaugedata

    clawdata.num_output_times = 30
    clawdata.tfinal = 3.0

    gaugedata.gauges = [
        [1, 1.0, 1.0, 0.0, 10.0],
        [2, 3.5, 0.5, 0.0, 10.0],
    ]

def configure_2d_forward(adjoint_output: Path):
    def _configure(runner):
        clawdata = runner.rundata.clawdata
        gaugedata = runner.rundata.gaugedata
        amrdata = runner.rundata.amrdata

        clawdata.num_output_times = 1
        clawdata.tfinal = 3.0

        # Test gauges
        gaugedata.gauges = []
        gaugedata.gauges.append([1, 1.0, 1.0, 0., 1e9])
        gaugedata.gauges.append([2, 3.5, 0.5, 0., 1e9])

        # AMR parameters
        amrdata.amr_levels_max = 2
        amrdata.refinement_ratios_x = [2]
        amrdata.refinement_ratios_y = [2]
        amrdata.refinement_ratios_t = [2]
        amrdata.flag_richardson_tol = 1e-5
        amrdata.flag2refine_tol = 0.02

        runner.rundata.adjointdata.adjoint_outdir = adjoint_output

    return _configure

@pytest.mark.regression
@pytest.mark.adjoint_forward
def test_acoustics_2d_adjoint_forward(tmp_path: Path, save: bool):
    example_path = Path(__file__).parent
    adjoint_path = example_path / "adjoint"
    adjoint_output = tmp_path / "_adjoint_output"

    run_example_for_test(
        test.AMRClawTestRunner,
        adjoint_output,
        adjoint_path,
        configure_runner=configure_2d_adjoint,
    )

    runner = run_example_for_test(
        test.AMRClawTestRunner,
        tmp_path,
        example_path,
        configure_runner=configure_2d_forward(adjoint_output),
    )

    runner.check_gauge(gauge_id=1, save=save)
    runner.check_gauge(gauge_id=2, save=save)

if __name__ == "__main__":
    pytest.main([__file__])
