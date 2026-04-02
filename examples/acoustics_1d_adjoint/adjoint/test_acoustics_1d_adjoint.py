#!/usr/bin/env python
"""
Test for the adjoint problem for 1D acoustics
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test

def set_adjoint_data(runner):
    """Set the adjoint data for the test case"""
    runner.set_data()

@pytest.mark.regression
@pytest.mark.adjoint
def test_acoustics_1d_adjoint(tmp_path: Path, save: bool):
    """Acoustics 1D adjoint test"""

    runner = test.AMRClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    
    runner.write_data()
    runner.build_executable()
    runner.run_code()

    runner.check_gauge(gauge_id=0)
    runner.check_gauge(gauge_id=1)


if __name__ == "__main__":
    pytest.main([__file__])
