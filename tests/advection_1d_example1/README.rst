
.. _amrclaw_tests_advection_1d_example1:

Advection 1D Example 1
------------------------------------------

1D advection with a constant velocity.  The equation is :math:`q_t + uq_x = 0`.


Quick regression test version:  regression_tests.py

* Runs only up to time t=0.2.
* Checks two gauges: Gauge 0 is at x=0.2 and Gauge 1 is at x=0.9.
