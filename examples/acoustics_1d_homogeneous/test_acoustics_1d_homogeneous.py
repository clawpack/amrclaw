#!/usr/bin/env python
"""
Regression tests for 1D acoustics in a homogeneous medium.
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test

def test_acoustics_1d_homogeneous(tmp_path: Path, save: bool):
    """Basic test for a 1D acoustics test case in a homogeneous medium"""

    runner = test.AMRClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    
    runner.set_data()
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 0.200000
    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gauges.append([0, 0.0, 0., 0.8])
    runner.rundata.gaugedata.gauges.append([1, 0.6, 0., 0.8])
    runner.write_data()

    runner.build_executable()

    runner.run_code()

    runner.check_gauge(save=save, gauge_id=0)
    runner.check_gauge(save=save, gauge_id=1)


if __name__ == "__main__":
    pytest.main([__file__])
