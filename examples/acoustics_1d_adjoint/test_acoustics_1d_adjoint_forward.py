"""
Regression tests for 1D acoustics with adjoint flagging.
"""

from pathlib import Path
import os
import pytest

import clawpack.amrclaw.test as test

# import adjoint.test_acoustics_1d_adjoint

def test_acoustics_1d_adjoint_forward(tmp_path: Path, save: bool):
    """Test for a 1D acoustics adjoint-flagging forward problem test case"""

    ctr = test.AMRClawTestRunner(tmp_path, test_path=__file__)
    
    # Run adjoint problem
    adjoint_output = ctr.temp_path / "_adjoint_output"
    
    if not adjoint_output.exists():
        os.makedirs(adjoint_output)
    adjoint_ctr = test.AMRClawTestRunner(adjoint_output, 
                                         test_path=ctr.test_path / "adjoint")
    adjoint_ctr.set_data(setrun_path=adjoint_ctr.test_path / "setrun.py")
    adjoint_ctr.write_data()
    adjoint_ctr.build_executable()
    adjoint_ctr.run_code()

    # Write problem data
    ctr.set_data()
    ctr.rundata.adjointdata.adjoint_outdir = adjoint_output
    ctr.write_data()

    ctr.build_executable()
    ctr.run_code()
    
    ctr.check_gauge(gauge_id=0)
    ctr.check_gauge(gauge_id=1)


if __name__ == "__main__":
    pytest.main([__file__])
