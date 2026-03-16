r"""
Defines the AMRClaw Clawpack Test Runner class for running PyTest based
regression tests in AMRClaw.

Refer to the documentation for PyTest to manage output and reporting.
"""

from pathlib import Path
import os
from typing import Optional

import clawpack.clawutil.test as test

# Clean library files whenever this module is used
if "CLAW" in os.environ:
    CLAW = Path(os.environ["CLAW"])
else:
    raise ValueError("Need to set CLAW environment variable.")

class AMRClawTestRunner(test.ClawpackTestRunner):

    def __init__(self, path: Path, test_path: Optional[Path]=None):
        super(AMRClawTestRunner, self).__init__(path, test_path=test_path)
        self.executable_name = 'xamr'


# Useful for running tests that need more than one example to be run, e.g. an
# adjoint test that needs to run the forward problem first.
def run_example_for_test(
    runner_cls,
    output_path: Path,
    test_path: Path,
    *,
    setrun_path: Path | None = None,
    configure_runner=None,
    build_kwargs: dict | None = None,
):
    """
    Build and run one example in a specified output directory.

    Parameters
    ----------
    runner_cls
        Runner class, e.g. AMRClawTestRunner.
    output_path
        Temporary directory for data files, executable, and output.
    test_path
        Path to the example directory.
    setrun_path
        Optional explicit path to setrun.py.
    configure
        Optional callback taking the runner after set_data() and before
        write_data().
    build_kwargs
        Optional keyword arguments passed to build_executable().
    """
    output_path.mkdir(parents=True, exist_ok=True)

    runner = runner_cls(output_path, test_path=test_path)
    runner.set_data(setrun_path=setrun_path)

    if configure_runner is not None:
        configure_runner(runner)

    runner.write_data()
    runner.build_executable(**(build_kwargs or {}))
    runner.run_code()
    return runner
