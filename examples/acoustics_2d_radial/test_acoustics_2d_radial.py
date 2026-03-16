#!/usr/bin/env python
"""
Regression tests for 2D acoustics.
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test

def test_acoustics_2d_radial(tmp_path: Path, save: bool):
    r"""Basic test for a 2D acoustics test case"""

    runner = test.AMRClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    
    runner.set_data()
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 0.2
    runner.rundata.regiondata.regions.append([1,1,0,1e10,-1.,1.,-1.,1])
    runner.rundata.regiondata.regions.append([1,3,0,1e10,-1.,1.,-.2,.2])
    runner.write_data()

    runner.build_executable()

    # Run code
    runner.run_code()

    # Perform tests
    runner.check_gauge(save=save, gauge_id=0)
    runner.check_gauge(save=save, gauge_id=1)
    runner.check_gauge(save=save, gauge_id=2)

if __name__ == "__main__":
    pytest.main([__file__])
