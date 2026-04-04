r"""
Defines the AMRClaw Clawpack Test Runner class for running PyTest based
regression tests in AMRClaw.

Refer to the documentation for PyTest to manage output and reporting.
"""

from pathlib import Path
import os
from typing import Optional

import clawpack.clawutil.test as test

# Set environment variable to avoid warning about missing CLAW variable
if "CLAW" in os.environ:
    CLAW = Path(os.environ["CLAW"])
else:
    raise ValueError("Need to set CLAW environment variable.")

class AMRClawTestRunner(test.ClawpackTestRunner):

    def __init__(self, path: Path, test_path: Optional[Path]=None):
        super(AMRClawTestRunner, self).__init__(path, test_path=test_path)
        self.executable_name = 'xamr'
