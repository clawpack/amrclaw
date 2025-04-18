"""
Regression tests for 1D acoustics in a homogeneous medium.
"""

import sys
import unittest

import clawpack.amrclaw.test as test


class Acoustics1DTest_Homogeneous(test.AMRClawRegressionTest):
    """Basic test for a 1D acoustics test case in a homogeneous medium"""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()
        
        self.rundata.clawdata.num_output_times = 1
        self.rundata.clawdata.tfinal = 0.200000

        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([0, 0.0, 0., 0.8])
        self.rundata.gaugedata.gauges.append([1, 0.6, 0., 0.8])

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=0)
        self.check_gauges(save=save, gauge_id=1)

        self.success = True



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Acoustics1DTest_Homogeneous()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()