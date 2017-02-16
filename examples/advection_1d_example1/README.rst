
.. advection_1d_example1:

Advection 1D Example 1
------------------------------------------

1D advection with a constant velocity.  The equation is :math:`q_t + uq_x = 0`.

The velocity is specified in setrun.py, along with the width parameter
of the Gaussian pressure pulse used for the initial condition.

Boundary conditions are periodic.

After running this code and creating plots via "make .plots", you
should be able to view the plots in _plots/_PlotIndex.html.

The true solution is computed in setplot.py.  If you change the initial
conditions, this should be changed accordingly.  (Except that if beta is 
modified and the code rerun, then then reloading setplot should read
in the latest value of beta used.)

Adapted from ‘$CLAW/classic/examples/advection_1d_example1’
to use adaptive mesh refinement.
