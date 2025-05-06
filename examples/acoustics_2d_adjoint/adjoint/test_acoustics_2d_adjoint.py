"""
Test for the adjoint problem for 2D acoustics
"""

import sys
import unittest

import clawpack.amrclaw.test as test
import clawpack.amrclaw.data as data

class Acoustics2DAdjointTest(test.AMRClawRegressionTest):

    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_output_times = 30
        self.rundata.clawdata.tfinal = 3.0

        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([1, 1.0, 1.0, 0., 10.])
        self.rundata.gaugedata.gauges.append([2, 3.5, 0.5, 0., 10.])

        self.write_rundata_objects()

        self.run_code()

        # Perform Tests
        self.check_gauges(save=save, gauge_id=1)
        self.check_gauges(save=save, gauge_id=2)

        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Acoustics2DAdjointTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()