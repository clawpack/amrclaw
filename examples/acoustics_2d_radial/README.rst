
.. _amrclaw_examples_acoustics_2d_radial:

Two-dimensional acoustics with radially symmetric initial data
==============================================================

Acoustics with radial symmetric initial conditions.  The solution should 
remain radially symmetric.  First run the code in the 1drad subdirectory to 
compute the "true solution" and then setplot.py contains code to produce a 
scatter plot of the computed 2d pressure vs. distance from origin and compare 
with the 1d solution.

Three gauges are specified, Gauge 0 is at the origin and Gauges 1 and 2 are
are both at radius 0.7, along the x-axis and along the diagonal.

