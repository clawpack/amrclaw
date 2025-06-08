"""
Test for the adjoint problem for 2D acoustics
"""

import sys
import unittest

import clawpack.amrclaw.test as test
import clawpack.amrclaw.data as data

class Acoustics1DAdjointTest(test.AMRClawRegressionTest):

    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()
        
        # self.rundata.clawdata.num_output_times = 10
        # self.rundata.clawdata.tfinal = 1.0

        # self.rundata.clawdata.use_fwaves = False

        self.write_rundata_objects()

        self.run_code()

        # Perform Tests
        self.check_gauges(save=save, gauge_id=0)
        self.check_gauges(save=save, gauge_id=1)

        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Acoustics1DAdjointTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()