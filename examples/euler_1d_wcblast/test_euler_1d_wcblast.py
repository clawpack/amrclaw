#!/usr/bin/env python
"""
Regression tests for 1D advection;.
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test

def test_euler_1d_wcblast(tmp_path: Path, save: bool):
    """Basic test for a 1D Euler case"""

    runner = test.AMRClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    
    runner.set_data()

    runner.rundata.clawdata.num_output_times = 2
    runner.rundata.clawdata.tfinal = 0.015

    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gauges.append([1, 0.3, 0, 1e9])
    runner.rundata.gaugedata.gauges.append([2, 0.85, 0, 1e9])
    runner.write_data()

    runner.build_executable()

    # Run code
    runner.run_code()

    # Perform tests
    runner.check_gauge(save=save, gauge_id=1)
    runner.check_gauge(save=save, gauge_id=2)

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
