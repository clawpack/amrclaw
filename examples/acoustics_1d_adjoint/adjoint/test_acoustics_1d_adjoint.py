"""
Test for the adjoint problem for 1D acoustics
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test

def test_acoustics_1d_adjoint(tmp_path: Path, save: bool):
    """Acoustics 1D adjoint test"""

    ctr = test.AMRClawTestRunner(tmp_path)
    ctr.set_data()
    ctr.write_data()
    ctr.build_executable()
    ctr.run_code()

    ctr.check_gauge(gauge_id=0)
    ctr.check_gauge(gauge_id=1)


if __name__ == "__main__":
    pytest.main([__file__])
