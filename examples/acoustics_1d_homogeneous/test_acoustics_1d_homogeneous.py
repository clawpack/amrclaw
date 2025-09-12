"""
Regression tests for 1D acoustics in a homogeneous medium.
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test

def test_acoustics_1d_homogeneous(tmp_path: Path, save: bool):
    """Basic test for a 1D acoustics test case in a homogeneous medium"""

    ctr = test.AMRClawTestRunner(tmp_path)
    
    ctr.set_data()
    ctr.rundata.clawdata.num_output_times = 1
    ctr.rundata.clawdata.tfinal = 0.200000
    ctr.rundata.gaugedata.gauges = []
    ctr.rundata.gaugedata.gauges.append([0, 0.0, 0., 0.8])
    ctr.rundata.gaugedata.gauges.append([1, 0.6, 0., 0.8])
    ctr.write_data()

    ctr.build_executable()

    ctr.run_code()

    ctr.check_gauge(gauge_id=0)
    ctr.check_gauge(gauge_id=1)


if __name__ == "__main__":
    pytest.main([__file__])
