"""
Regression tests for swril 3D advection.
"""

from __future__ import absolute_import
import sys
import os
import unittest

import clawpack.amrclaw.test as test


class Advection3DSwirlTest(test.AMRClawRegressionTest):
    r"""Basic test for a 3D advection test case"""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_output_times = 2
        self.rundata.clawdata.tfinal = .1

        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([1, 0.55, 0.4, 0.4, 0., 1e9])
        self.rundata.gaugedata.gauges.append([2, 0.45, 0.6, 0.4, 0., 1e9])

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_frame(save=save, frame_num=1, file_name="regression_data_test2.txt")
        self.check_frame(save=save, frame_num=2, file_name="regression_data_test3.txt")

        self.check_gauges(save=save, gauge_id=1)
        self.check_gauges(save=save, gauge_id=2)

        self.success = True



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Advection3DSwirlTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()