"""
Regression tests for swril 3D advection.
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test

def test_advection_3d_swirl(tmp_path: Path, save: bool):
    r"""Basic test for a 3D advection on a swirl test case"""

    runner = test.AMRClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    
    runner.set_data()
    runner.rundata.clawdata.num_output_times = 2
    runner.rundata.clawdata.tfinal = .1
    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gauges.append([1, 0.55, 0.4, 0.4, 0., 1e9])
    runner.rundata.gaugedata.gauges.append([2, 0.45, 0.6, 0.4, 0., 1e9])
    runner.write_data()

    runner.build_executable()

    runner.run_code()

    runner.check_gauge(gauge_id=1)
    runner.check_gauge(gauge_id=2)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
