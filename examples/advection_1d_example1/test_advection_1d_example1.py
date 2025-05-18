"""
Regression tests for 1D advection;.
"""

import sys
import unittest

import clawpack.amrclaw.test as test


class Advection1DTest(test.AMRClawRegressionTest):
    """Basic test for a 1D acoustics advection case"""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_cells[0] = 40
        self.rundata.clawdata.num_output_times = 1
        self.rundata.clawdata.tfinal = 0.200000

        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([0, 0.2, 0, 1e9])
        self.rundata.gaugedata.gauges.append([1, 0.9, 0, 1e9])

        self.rundata.amrdata.refinement_ratios_x = [2, 2]
        self.rundata.amrdata.refinement_ratios_t = [2, 2]
        self.rundata.amrdata.flag2refine_tol = 0.1
        self.rundata.amrdata.regrid_buffer_width = 2

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
            test = Advection1DTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()