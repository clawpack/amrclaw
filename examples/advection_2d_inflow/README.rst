
.. _amrclaw_examples_advection_2d_inflow:

Two-dimensional advection with inflow boundaries
================================================

See `bc2amr.f90` for the modified routine that implements inflow boundary
conditions.  This routine and `qinit.f` use `qtrue.f90` to specify the true
solution for initial and boundary conditions.

The advection velocities `u` and `v` can be changed in `setrun.py`, but they
are assumed to both be positive since the inflow boundary conditions are
implemented only at the left and bottom boundaries.
