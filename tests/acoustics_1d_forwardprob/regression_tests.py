"""
Regression tests for 1D acoustics with adjoint flagging.
"""

from __future__ import print_function
from __future__ import absolute_import
import sys,os
import unittest

thisfile = os.path.realpath(__file__)
testdir = os.path.split(thisfile)[0]

import clawpack.amrclaw.test as test


class Acoustics1DAdjointTest(test.AMRClawRegressionTest):
    r"""Basic test for a 1D acoustics adjoint-flagging forward problem test case"""


    def runTest(self, save=False):
        
        # Run adjoint problem
        adjointdir = testdir + '/adjoint'

        # Running the adjoint problem
        os.chdir(adjointdir)
        os.system('make -s new')
        os.system('make .output > /dev/null')
        os.chdir(testdir)

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=0)
        self.check_gauges(save=save, gauge_id=1)
        self.check_gauges(save=save, gauge_id=2)

        # self.check_gauges(save=save, gauge_num=1,
        #                   regression_data_path='regression_data_test2.txt')
        # self.check_gauges(save=save, gauge_num=2,
        #                   regression_data_path='regression_data_test3.txt')

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
