
"""
Run several regression tests with different parameters.
"""

from setrun_regression import setrun
from setplot import setplot
from clawpack.clawutil.runclaw import runclaw
from clawpack.visclaw.plotclaw import plotclaw
from clawpack.clawutil.compare_regression_tests import compare_regression_tests
import os,sys

if 1:
    # Compile code from scratch to make sure up to date:
    os.system('make new')

compare_results = True  # will attempt to download archived results to compare

regression_dir = "_regression_tests"

if os.path.exists(regression_dir):
    ans = raw_input("Directory %s exists, ok to overwrite? " % regression_dir)
    if ans=='y':
        os.system('rm -rf %s' % regression_dir)
    else:
        print "*** Aborting regression tests"
        sys.exit()

os.system('mkdir %s' % regression_dir)

# initialize rundata using setrun but then change some things for each run:
rundata = setrun()
clawdata = rundata.clawdata
amrdata = rundata.amrdata


#--------------
# Test 1:
#--------------

suffix = "_test1"
outdir = regression_dir + "/_output" + suffix
plotdir = regression_dir + "/_plots" + suffix

amrdata.amr_levels_max = 1
rundata.write()

runclaw(xclawcmd = "xamr", outdir=outdir, print_git_status=True)
plotclaw(outdir=outdir, plotdir=plotdir)

#--------------
# Test 2:
#--------------
amrdata.amr_levels_max = 2
rundata.write()

suffix = "_test2"
outdir = regression_dir + "/_output" + suffix
plotdir = regression_dir + "/_plots" + suffix

runclaw(xclawcmd = "xamr", outdir=outdir, print_git_status=True)
plotclaw(outdir=outdir, plotdir=plotdir)

#--------------
# Test 3:
#--------------

amrdata.amr_levels_max = 3
rundata.write()

suffix = "_test3"
outdir = regression_dir + "/_output" + suffix
plotdir = regression_dir + "/_plots" + suffix

runclaw(xclawcmd = "xamr", outdir=outdir, print_git_status=True)
plotclaw(outdir=outdir, plotdir=plotdir)

print "Output and plots are in ", regression_dir

if compare_results:
    compare_regression_tests(regression_dir)

