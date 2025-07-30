"""
Regression tests for 2D advection on a square.
"""

import sys
import unittest

import clawpack.amrclaw.test as test
import clawpack.amrclaw.data as data


class Advection2DSquareTest(test.AMRClawRegressionTest):
    r"""Basic test for a 2D advection test case"""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()

        self.rundata.clawdata.num_output_times = 1
        self.rundata.clawdata.tfinal = 0.4

        self.rundata.gaugedata.gauges = []
        self.rundata.gaugedata.gauges.append([1, 0.65, 0.4, 0., 10.])
        self.rundata.gaugedata.gauges.append([2, 0.2, 0.8, 0., 10.])

        # Test newer flagregions

        # The entire domain restricted to level 1 for illustration:
        # Note that this is a rectangle specified in the new way:
        # (other regions below will force/allow more refinement)
        flagregion = data.FlagRegion(num_dim=2)
        flagregion.name = 'Region_domain'
        flagregion.minlevel = 1
        flagregion.maxlevel = 1
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [0.,1.,0.,1.]  # = [x1,x2,y1,y2]
        self.rundata.flagregiondata.flagregions.append(flagregion)

        # Another rectangle specified in the new way:
        flagregion = data.FlagRegion(num_dim=2)
        flagregion.name = 'Region_3levels'
        flagregion.minlevel = 1
        flagregion.maxlevel = 3
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [0.,1.,0.,0.7]
        self.rundata.flagregiondata.flagregions.append(flagregion)

        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
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
            test = Advection2DSquareTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()