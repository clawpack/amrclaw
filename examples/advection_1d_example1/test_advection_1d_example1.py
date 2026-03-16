#!/usr/bin/env python
"""
Regression tests for 1D advection;.
"""

<<<<<<< HEAD
import sys
import unittest
import numpy
=======
from pathlib import Path
import pytest
>>>>>>> 83d5ea6 (Convert over existing unittest tests)

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

    runner.rundata.amrdata.refinement_ratios_x = [2, 2]
    runner.rundata.amrdata.refinement_ratios_t = [2, 2]
    runner.rundata.amrdata.flag2refine_tol = 0.1
    runner.rundata.amrdata.regrid_buffer_width = 2
    runner.write_data()

    runner.build_executable()

    runner.run_code()

<<<<<<< HEAD
        # Set up 3 gauges at the same location with different formats
        gauge_location = 0.2
        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([0, gauge_location, 0, 1e9])  # ascii
        self.rundata.gaugedata.gauges.append([1, gauge_location, 0, 1e9])  # binary32
        self.rundata.gaugedata.gauges.append([2, gauge_location, 0, 1e9])  # binary64

        # Set per-gauge formats: ascii, binary32, binary64
        self.rundata.gaugedata.file_format = ['ascii', 'binary32', 'binary64']

        self.rundata.amrdata.refinement_ratios_x = [2, 2]
        self.rundata.amrdata.refinement_ratios_t = [2, 2]
        self.rundata.amrdata.flag2refine_tol = 0.1
        self.rundata.amrdata.regrid_buffer_width = 2

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Load all three gauge files
        gauge_ascii = gauges.GaugeSolution(0, path=self.temp_path)
        gauge_binary32 = gauges.GaugeSolution(1, path=self.temp_path)
        gauge_binary64 = gauges.GaugeSolution(2, path=self.temp_path)

        # Compare all formats against each other
        # Use ASCII as reference for comparisons
        indices = [0]  # Compare first component (q[0])
        
        # binary32 vs ascii: binary32 has ~7 decimal digits precision
        numpy.testing.assert_allclose(
            gauge_binary32.q[indices[0], :],
            gauge_ascii.q[indices[0], :],
            rtol=1e-7, atol=1e-7,
            err_msg="binary32 gauge does not match ascii gauge")
        
        # binary64 vs ascii: binary64 has ~15-17 decimal digits precision
        numpy.testing.assert_allclose(
            gauge_binary64.q[indices[0], :],
            gauge_ascii.q[indices[0], :],
            rtol=1e-15, atol=1e-15,
            err_msg="binary64 gauge does not match ascii gauge")
        
        # binary32 vs binary64: limited by binary32 precision
        numpy.testing.assert_allclose(
            gauge_binary32.q[indices[0], :],
            gauge_binary64.q[indices[0], :],
            rtol=1e-7, atol=1e-7,
            err_msg="binary32 gauge does not match binary64 gauge")

        self.success = True



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Advection1DTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
=======
    runner.check_gauge(save=save, gauge_id=0)
    runner.check_gauge(save=save, gauge_id=1)

if __name__ == "__main__":
    pytest.main([__file__])
>>>>>>> 83d5ea6 (Convert over existing unittest tests)
