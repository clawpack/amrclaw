This directory contains a few examples to get users started using
this software.  In most directories you should be able to type:
    make .plots
and the code should run and a set of plots produced.  However, see the
individual README files in the directories.

The examples in this directory can also be used as regression tests.
Each example is run and the set of plots are then compared with those
archived in the gallery that can be found online in the documentation.
To do this comparison, you must first clone the repository 
  https://github.com/clawpack/clawpack.github.com

To run all tests:
    python run_tests.py
or: 
    make tests

To clean up afterwards, removing all executables, output, and test results:
    make clobber

It may take some time for these examples to all run.

A quicker set of tests can be run from the amrclaw/tests directory.  See the
README file in that directory.

