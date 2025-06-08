"""
Regression tests for 2D acoustics.
"""

import sys
import unittest

import clawpack.amrclaw.test as test


class Acoustics2DRadialTest(test.AMRClawRegressionTest):
    r"""Basic test for a 2D acoustics test case"""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_output_times = 1
        self.rundata.clawdata.tfinal = 0.2
        self.rundata.regiondata.regions.append([1,1,0,1e10,-1.,1.,-1.,1])
        self.rundata.regiondata.regions.append([1,3,0,1e10,-1.,1.,-.2,.2])

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=0)
        self.check_gauges(save=save, gauge_id=1)
        self.check_gauges(save=save, gauge_id=2)

        self.success = True



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Acoustics2DRadialTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()