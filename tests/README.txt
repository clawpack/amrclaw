This directory is for regression tests.

To run all tests:
    python run_tests.py
or: 
    make tests

To clean up afterwards, removing all executables, output, and test results:
    make clobber

Each test does a short run with regions and gauges set to exercise the code,
The test passes if the code runs and if the sum of t values and of q values
agree with archived results for each gauge.


Developers: To create new archived results for a test case, modify 
the main program in the file regression_tests.py for this test to set 
    save_new_regression_data = True
and then execute as a script.  Then 'git add' and issue a pull request if
you believe the new results are more correct.
