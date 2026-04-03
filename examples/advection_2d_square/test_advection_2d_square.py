#!/usr/bin/env python
"""
Regression tests for 2D advection on a square.
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test
import clawpack.amrclaw.data as data

def test_advection_2d_square(tmp_path: Path, save: bool):
    """Basic test for a 2D advection test case"""

    runner = test.AMRClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    
    runner.set_data()
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 0.4
    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gauges.append([1, 0.65, 0.4, 0., 10.])
    runner.rundata.gaugedata.gauges.append([2, 0.2, 0.8, 0., 10.])

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
    runner.rundata.flagregiondata.flagregions.append(flagregion)

    # Another rectangle specified in the new way:
    flagregion = data.FlagRegion(num_dim=2)
    flagregion.name = 'Region_3levels'
    flagregion.minlevel = 1
    flagregion.maxlevel = 3
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [0.,1.,0.,0.7]
    runner.rundata.flagregiondata.flagregions.append(flagregion)

    runner.write_data()

    runner.build_executable()

    runner.run_code()

    runner.check_gauge(save=save, gauge_id=1)
    runner.check_gauge(save=save, gauge_id=2)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
