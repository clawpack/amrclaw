#!/usr/bin/env python
"""
Regression tests for 1D acoustics with adjoint flagging.
"""

from pathlib import Path
import pytest

import clawpack.amrclaw.test as test

def test_acoustics_1d_adjoint_forward(tmp_path: Path, save: bool):
    """Test for a 1D acoustics adjoint-flagging forward problem test case"""
    
    example_path = Path(__file__).parent
    adjoint_path = example_path / "adjoint"
    adjoint_output = tmp_path / "_adjoint_output"

    test.run_example_for_test(
        test.AMRClawTestRunner,
        adjoint_output,
        adjoint_path,
    )

    runner = test.run_example_for_test(
        test.AMRClawTestRunner,
        tmp_path,
        example_path,
        configure_runner=lambda r: setattr(r.rundata.adjointdata,
                                           "adjoint_outdir", adjoint_output),
    )

    runner.check_gauge(0, save=save)
    runner.check_gauge(1, save=save)

if __name__ == "__main__":
    pytest.main([__file__])
