
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

Extrapolation BCs
------------------
The code is set up to use extrapolation boundary conditions at all
boundaries. This does a reasonably good job of providing non-reflecting 
boundaries, but there are some artifacts visible at later times.

Absorbing boundary layer
------------------------
New cababilities have been added to Clawpack 5.5.0 to allow extending the
computational domain with an aborbing boundary layer that does a better job
of eliminating artificial reflections.  [Add more discussion.]
To try this version::
    make all -f Makefile_abl

