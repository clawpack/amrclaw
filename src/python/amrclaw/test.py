r"""
Defines the Classic Clawpack Test Runner class for running PyTest based
regression tests in classic clawpack.

Refer to the documentation for PyTest to manage output and reporting.
"""

from pathlib import Path
import os

import clawpack.clawutil.test as test

# Clean library files whenever this module is used
if "CLAW" in os.environ:
    CLAW = Path(os.environ["CLAW"])
else:
    raise ValueError("Need to set CLAW environment variable.")

for lib_path in (CLAW / "amrclaw" / "src" / "1d").glob("*.o"):
    lib_path.unlink()
for lib_path in (CLAW / "amrclaw" / "src" / "2d").glob("*.o"):
    lib_path.unlink()
for lib_path in (CLAW / "amrclaw" / "src" / "3d").glob("*.o"):
    lib_path.unlink()

class AMRClawTestRunner(test.ClawpackTestRunner):

    def __init__(self, path: Path):
        super(AMRClawTestRunner, self).__init__(path)
        self.executable_name = 'xamr'
