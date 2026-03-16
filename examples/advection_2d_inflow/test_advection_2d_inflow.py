#!/usr/bin/env python
"""
Regression tests for 2D advection on an annulus.
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test

def test_advection_2d_inflow(tmp_path: Path, save: bool):
    r"""Basic test for a 2D advection on an annulus test case"""

    runner = test.AMRClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    
    runner.set_data()
    runner.rundata.probdata.add_param('u',     0.5,  'ubar advection velocity')
    runner.rundata.probdata.add_param('v',     1.0,  'vbar advection velocity')
    runner.rundata.clawdata.output_style = 3
    runner.rundata.clawdata.output_step_interval = 10
    runner.rundata.clawdata.total_steps = 10
    runner.rundata.clawdata.bc_lower[1] = 'extrap' # opposed to user
    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gauges.append([1, 0.08, 0.9, 0., 10.])
    runner.rundata.gaugedata.gauges.append([2, 0.06, 0.391, 0., 10.])
    runner.write_data()

    runner.build_executable(make_vars={'TEST': True}, verbose=True)

    runner.run_code()

    runner.check_gauge(save=save, gauge_id=1)
    runner.check_gauge(save=save, gauge_id=2)

if __name__ == "__main__":
    pytest.main([__file__])
