
.. _amrclaw_tests_advection_2d_square:

Two-dimensional advection of a square pulse 
===========================================

With periodic boundary conditions.

Quick regression test version:  regression_tests.py

* Allows 3 levels only for Y <= 0.7, elsewhere 1 level.
* Runs only up to time t=0.4.
* Checks sum of t values and q values at two different gauges.
  Tests pass if these agree with specified values in regression_tests.py.
