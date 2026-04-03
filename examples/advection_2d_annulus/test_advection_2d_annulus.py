#!/usr/bin/env python
"""
Regression tests for 2D advection on an annulus.
"""

from pathlib import Path
import pytest

import numpy as np

import clawpack.amrclaw.test as test

def test_advection_2d_annulus(tmp_path: Path, save: bool):
    r"""Basic test for a 2D advection on an annulus test case"""

    runner = test.AMRClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    
    runner.set_data()
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 0.500000
    runner.rundata.regiondata.regions.append([1, 2, 0.0, 10.0, 0.2, 1.0, 0.0, 2.0 * np.pi])
    runner.rundata.regiondata.regions.append([3, 3, 0.0, 10.0, 0.5, 1.0, 0.0, 0.5 * np.pi])

    runner.write_data()

    runner.build_executable()

    runner.run_code()

    runner.check_gauge(save=save, gauge_id=1)
    runner.check_gauge(save=save, gauge_id=2)

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
