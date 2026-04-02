#!/usr/bin/env python
"""
Test for the adjoint problem for 2D acoustics
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

@pytest.mark.regression
@pytest.mark.adjoint
def test_acoustics_2d_adjoint(tmp_path: Path, save: bool):
    adjoint_path = Path(__file__).parent

    runner = run_example_for_test(
        test.AMRClawTestRunner,
        tmp_path,
        adjoint_path,
        configure_runner=configure_2d_adjoint,
    )

    runner.check_gauge(gauge_id=1, save=save)
    runner.check_gauge(gauge_id=2, save=save)

if __name__ == "__main__":
    pytest.main([__file__])
