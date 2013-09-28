

.. _amrclaw_tests_advection_2d_annulus:

Two-dimension advection around an annulus, on a mapped grid
===========================================================

Quick regression test version:  regression_tests.py

* Forces 3 levels for 0 <= theta <= pi/2, elsewhere allows 2 levels.
* Runs only up to time t=1.0.
* Checks sum of t values and q values at two different gauges.
  Tests pass if these agree with specified values in regression_tests.py.

