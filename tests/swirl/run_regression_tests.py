
"""
Run several regression tests with different parameters.

Check the settings in the __main__ program to specify what is done.

To modify for a different app or set of tests, see the code labelled "Test 1", 
"Test 2", etc.
"""

from setrun_regression import setrun
from setplot import setplot
from clawpack.clawutil.runclaw import runclaw
from clawpack.visclaw.plotclaw import plotclaw
from clawpack.clawutil.compare_regression_tests import compare_regression_tests
import os,sys


def run_regression_tests(regression_dir="_regression_tests", \
            regression_output_files="all", \
            regression_plot_files="all", \
            make_new=True, run_tests=True, compare_results=True, \
            relocatable=False):

    if make_new:
        # Compile code from scratch to make sure up to date:
        os.system('make new')

    if run_tests:

        if os.path.exists(regression_dir):
            ans = raw_input("Directory %s exists, ok to overwrite? " % regression_dir)
            if ans=='y':
                os.system('rm -rf %s' % regression_dir)
            else:
                print "*** Aborting regression tests"
                sys.exit()

        os.system('mkdir %s' % regression_dir)

        # -----------------------------------------------------------
        # Define the regression test runs:
        # -----------------------------------------------------------

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

        # ----------------------------------------------------
        # End of test case definitions
        # ----------------------------------------------------

        print "Output and plots are in ", regression_dir

    regression_ok = None

    if compare_results:
        regression_ok = compare_regression_tests(regression_dir,\
                        regression_output_files, regression_plot_files, \
                        relocatable=relocatable)
        if regression_ok:
            print "The specified regression files are identical in all"
            print "    directories checked"
        else:
            print "*** Some regression tests did not pass"

    return regression_ok

if __name__=="__main__":

    # Location for regression results:
    regression_dir = "_regression_tests"

    # Which output files to check in order to "pass" this test:
    # If these are identical for each test then run_regression_tests returns True.
    regression_output_files = ['fort.q0002','fort.gauge']
    regression_plot_files = ['frame0002fig0.png','gauge0001fig300.png']

    # Specify what to do:
    make_new = False        # run 'make new' to recompile everything?
    run_tests = False       # run all the tests?
    compare_results = True  # download archived results and compare?
    relocatable = True      # copy all image files so regression_dir can
                            #      be moved, e.g. posted for discussion

    run_regression_tests(regression_dir, regression_output_files, \
                regression_plot_files, \
                make_new, run_tests, compare_results, \
                relocatable=relocatable)
