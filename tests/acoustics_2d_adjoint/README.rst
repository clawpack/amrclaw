
.. _amrclaw_tests_acoustics_2d_adjoint:

Two-dimensional acoustics adjoint-flagging forward problem 
===========================================================

Quick regression test version:  regression_tests.py

* Runs only up to time t=0.5.
* Checks sum of t values and q values at two different gauges.
  Tests pass if these agree with specified values in regression_tests.py.
* Only tests the forward problem.  The adjoint problem is first solved
  by running `make .output` in the adjoint subdirectory, but no tests are
  performed on that output (although adjoint/regression_data does contain
  the expected gauge output for comparison).
