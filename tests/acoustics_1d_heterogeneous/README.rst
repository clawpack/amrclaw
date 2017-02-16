
.. _amrclaw_tests_acoustics_1d_heterogeneous:

Acoustics 1D -- Heterogeneous example
------------------------------------------

Two gauges are specified, Gauge 0 is at the origin and Gauge 1 is
at -2, along the x-axis.

Quick regression test version:  regression_tests.py

* Runs only up to time t=0.2.
* Checks sum of t values and q values at three different gauges.
  Tests pass if these agree with specified values in regression_tests.py.

