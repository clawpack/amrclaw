#!/usr/bin/env python
"""
Regression tests for 2D advection on an annulus.
"""

import sys
import os
import unittest
import subprocess
import inspect
import shutil

import clawpack.amrclaw.test as test


class Advection2DBoundaryTest(test.AMRClawRegressionTest):
    r"""Basic test for a 2D advection on an annulus test case"""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()

        self.rundata.probdata.add_param('u',     0.5,  'ubar advection velocity')
        self.rundata.probdata.add_param('v',     1.0,  'vbar advection velocity')

        self.rundata.clawdata.output_style = 3
        self.rundata.clawdata.output_step_interval = 10
        self.rundata.clawdata.total_steps = 10

        self.rundata.clawdata.bc_lower[1] = 'extrap' # opposed to user

        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([1, 0.08, 0.9, 0., 10.])
        self.rundata.gaugedata.gauges.append([2, 0.06, 0.391, 0., 10.])

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1)
        self.check_gauges(save=save, gauge_id=2)

        self.success = True


    def build_executable(self, executable_name="xamr"):
        r"""Build executable by running `make .exe` in test directory.

        Moves the resulting executable to the temporary directory.


        """

        try:
            self.stdout.write("Test path and class info:\n")
            self.stdout.write("  class: %s\n" % str(self.__class__))
            self.stdout.write("  class file: %s\n" % str(inspect.getfile(self.__class__)))
            self.stdout.write("  test path: %s\n" % str(self.test_path))
            self.stdout.write("  temp path: %s\n" % str(self.temp_path))
            subprocess.check_call("cd %s ; make TEST=T .exe" % self.test_path, 
                                                               stdout=self.stdout,
                                                               stderr=self.stderr,
                                                               shell=True)
        except subprocess.CalledProcessError as e:
            self.tearDown()
            raise e

        self.executable_name = executable_name
        shutil.move(os.path.join(self.test_path, self.executable_name),  
                    self.temp_path)


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Advection2DBoundaryTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()