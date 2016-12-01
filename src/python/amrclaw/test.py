r"""
Execute nosetests in all subdirectories, to run a series of quick
regression tests.

Sends output and result/errors to separate files to simplify checking
results and looking for errors.
"""

from __future__ import absolute_import
import os
import glob

import clawpack.clawutil.test
import clawpack.pyclaw.util

# Clean library files whenever this module is used
if "CLAW" in os.environ:
    CLAW = os.environ["CLAW"]
else:
    raise ValueError("Need to set CLAW environment variable.")

for lib_path in [os.path.join(CLAW,"amrclaw","src","2d"),
                 os.path.join(CLAW,"amrclaw","src","3d")]:
    for path in glob.glob(os.path.join(lib_path,"*.o")):
        os.remove(path)
    for path in glob.glob(os.path.join(lib_path,"*.mod")):
        os.remove(path)


class AMRClawRegressionTest(clawpack.clawutil.test.ClawpackRegressionTest):

    r"""Base AMRClaw regression test setup derived from ClawpackRegressionTest

    """

    __doc__ += clawpack.pyclaw.util.add_parent_doc(
                                  clawpack.clawutil.test.ClawpackRegressionTest)


    def build_executable(self, executable_name="xamr"):
        r"""Build executable by running `make .exe` in test directory.

        Moves the resulting executable to the temporary directory.


        """

        super(AMRClawRegressionTest, self).build_executable(
                                                executable_name=executable_name)
