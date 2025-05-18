"""
Regression tests for 2D acoustics with adjoint flagging.
"""

from pathlib import Path
import sys
import shutil
import unittest

import clawpack.amrclaw.test as test

from adjoint.test_acoustics_2d_adjoint import Acoustics2DAdjointTest 

class Acoustics2DAdjointForwardTest(test.AMRClawRegressionTest):
    r"""Basic test for a 2D acoustics adjoint-flagging forward problem test case"""


    def runTest(self, save=False):
        
        # Run adjoint problem
        try:
            adjoint_run = Acoustics2DAdjointTest()    
            adjoint_run.setUp()
            adjoint_run.runTest()
            
            # Copy output to local directory
            adjoint_output = Path(self.temp_path) / "_adjoint_output"

            if Path(adjoint_output).exists():
                shutil.rmtree(adjoint_output)
            shutil.copytree(adjoint_run.temp_path, adjoint_output)
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

        # Look for adjoint data
        self.rundata.adjointdata.adjoint_outdir = adjoint_output.resolve()

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
            test = Acoustics2DAdjointForwardTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
