"""
Regression tests for 2D acoustics with adjoint flagging.
"""

import sys
# import os
import unittest

import clawpack.amrclaw.test as test
import clawpack.clawutil.runclaw

adjoint_path = os.path.abspath(os.path.join(self.test_path, "adjoint"))

class Acoustics2DAdjointRun(test.AMRClawRegressionTest):

    def __init__(self, methodName="runTest"):

        super(test.AMRClawRegressionTest, self).__init__(methodName=methodName)

        self.test_path = adjoint_path

    def runTest(self, save=False):

        
        # Write out data files
        self.load_rundata()

        adjoint_setrun.rundata.clawdata.num_output_times = 30
        adjoint_setrun.rundata.clawdata.tfinal = 3.0

        adjoint_setrun.rundata.gaugedata.gauges.append([1, 1.0, 1.0, 0., 10.])
        adjoint_setrun.rundata.gaugedata.gauges.append([2, 3.5, 0.5, 0., 10.])

        self.write_rundata_objects()

        

class Acoustics2DAdjointTest(test.AMRClawRegressionTest):
    r"""Basic test for a 2D acoustics adjoint-flagging forward problem test case"""


    def runTest(self, save=False):
        
        # Run adjoint problem
        adjoint_run = Acoustics2DAdjointRun()    
        try:
            adjoint_run.setUp()
            adjoint_run.runTest(save=True)
        finally:
            adjoint_run.tearDown()

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_output_times = 1
        self.rundata.clawdata.tfinal = 3.0

        # Test gauges
        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([1, 1.0, 1.0, 0., 1e9])
        self.rundata.gaugedata.gauges.append([2, 3.5, 0.5, 0., 1e9])

        # AMR parameters
        self.rundata.amrdata.amr_levels_max = 2
        self.rundata.amrdata.refinement_ratios_x = [2]
        self.rundata.amrdata.refinement_ratios_y = [2]
        self.rundata.amrdata.refinement_ratios_t = [2]
        self.rundata.amrdata.flag_richardson_tol = 1e-5
        self.rundata.amrdata.flag2refine_tol = 0.02

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
            test = Acoustics2DAdjointTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
