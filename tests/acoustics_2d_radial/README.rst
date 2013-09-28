
.. _amrclaw_tests_acoustics_2d_radial:

Two-dimensional acoustics with radially symmetric initial data
==============================================================

Three gauges are specified, Gauge 0 is at the origin and Gauges 1 and 2 are
are both at radius 0.7, along the x-axis and along the diagonal.

Quick regression test version:  regression_tests.py

* Allows 3 levels for -0.2 <= y <= 0.2, elsewhere allows 1 level.
* Runs only up to time t=0.2.
* Checks sum of t values and q values at three different gauges.
  Tests pass if these agree with specified values in regression_tests.py.

