"""
Regression tests for 1D advection;.
"""

import sys
import unittest

import clawpack.amrclaw.test as test


class Euler1DTest(test.AMRClawRegressionTest):
    """Basic test for a 1D Euler case"""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_output_times = 2
        self.rundata.clawdata.tfinal = 0.015

        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([1, 0.3, 0, 1e9])
        self.rundata.gaugedata.gauges.append([2, 0.85, 0, 1e9])

        # amrdata.max1d = 500

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
            test = Euler1DTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
