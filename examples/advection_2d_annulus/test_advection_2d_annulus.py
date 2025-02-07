"""
Regression tests for 2D advection on an annulus.
"""

import sys
import unittest

import numpy as np

import clawpack.amrclaw.test as test


class Advection2DAnnulusTest(test.AMRClawRegressionTest):
    r"""Basic test for a 2D advection on an annulus test case"""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_output_times = 1
        self.rundata.clawdata.tfinal = 0.500000

        self.rundata.regiondata.regions.append([1, 2, 0.0, 10.0, 0.2, 1.0, 0.0, 2.0 * np.pi])
        self.rundata.regiondata.regions.append([3, 3, 0.0, 10.0, 0.5, 1.0, 0.0, 0.5 * np.pi])

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1)
        self.check_gauges(save=save, gauge_id=2)

        self.success = True



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Advection2DAnnulusTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()