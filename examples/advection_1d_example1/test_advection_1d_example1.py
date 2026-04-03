#!/usr/bin/env python
"""
Regression tests for 1D advection;.
"""

from pathlib import Path
import pytest
import numpy as np

import clawpack.amrclaw.test as test
import clawpack.pyclaw.gauges as gauges

def test_advection_1d_example1(tmp_path: Path, save: bool):
    r"""Basic test for a 1D advection test case"""

    runner = test.AMRClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    
    runner.set_data()
    runner.rundata.clawdata.num_cells[0] = 40
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 0.200000

    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gauges.append([0, 0.2, 0, 1e9])
    runner.rundata.gaugedata.gauges.append([1, 0.9, 0, 1e9])
    runner.rundata.gaugedata.gauges.append([2, 0.2, 0, 1e9])
    runner.rundata.gaugedata.gauges.append([3, 0.2, 0, 1e9])
    
    runner.rundata.gaugedata.file_format = ['ascii', 'ascii', 
                                            'binary32', 'binary64']

    runner.rundata.amrdata.refinement_ratios_x = [2, 2]
    runner.rundata.amrdata.refinement_ratios_t = [2, 2]
    runner.rundata.amrdata.flag2refine_tol = 0.1
    runner.rundata.amrdata.regrid_buffer_width = 2
    runner.write_data()

    runner.build_executable()

    runner.run_code()

    # Original gauge check
    runner.check_gauge(save=save, gauge_id=1)

    # Testing for gauge formatting: load all three gauge files
    gauge_ascii = gauges.GaugeSolution(0, path=runner.temp_path)
    gauge_binary32 = gauges.GaugeSolution(2, path=runner.temp_path)
    gauge_binary64 = gauges.GaugeSolution(3, path=runner.temp_path)

    # Compare all formats against each other
    # Use ASCII as reference for comparisons
    indices = [0]  # Compare first component (q[0])
    
    # binary32 vs ascii: binary32 has ~7 decimal digits precision
    np.testing.assert_allclose(
        gauge_binary32.q[indices[0], :],
        gauge_ascii.q[indices[0], :],
        rtol=1e-7, atol=1e-7,
        err_msg="binary32 gauge does not match ascii gauge")
    
    # binary64 vs ascii: binary64 has ~15-17 decimal digits precision
    np.testing.assert_allclose(
        gauge_binary64.q[indices[0], :],
        gauge_ascii.q[indices[0], :],
        rtol=1e-15, atol=1e-15,
        err_msg="binary64 gauge does not match ascii gauge")
    
    # binary32 vs binary64: limited by binary32 precision
    np.testing.assert_allclose(
        gauge_binary32.q[indices[0], :],
        gauge_binary64.q[indices[0], :],
        rtol=1e-7, atol=1e-7,
        err_msg="binary32 gauge does not match binary64 gauge")

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
