
.. _amrclaw_tests_advection_1d_example1:

Advection 1D Example 1
------------------------------------------

1D advection with a constant velocity.  The equation is :math:`q_t + uq_x = 0`.

Two gauges are specified, Gauge 0 is at the origin and Gauge 1 is
at 1, along the x-axis.

Quick regression test version:  regression_tests.py

* Runs only up to time t=0.2.
* Checks sum of t values and q values at three different gauges.
  Tests pass if these agree with specified values in regression_tests.py.
