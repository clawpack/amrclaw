
.. euler_1d_wcblast:

Euler equations in 1d -- Woodward-Colella blast wave problem
------------------------------------------------------------

The Woodward-Colella interacting blast wave problem consists of the
one-dimensional Euler equations of compressible gas dynamics together with
initial data containing two discontinuities in pressure that lead to strong
shock waves that interact with one another.

This code is based on the discussion in Example 15.1 of 
`Finite Volume Methods for Hyperbolic Problems
<http://www.clawpack.org/book.html>`_ by R. J. LeVeque.
Figures 15.5 and 15.6 show the solution at time t = 0.038 as computed with
various methods.

To better illustrate the AMR solution, in this example the solution on each
grid is plotted as a piecewise constant function over each grid cell, rather
than as a piecewise linear function connecting cell averages (as is usually
done).  This is specified in `setplot.py` by using the new `ClawPlotItem`
option `plot_type = 1d_pwconst`, introduced in Version 5.4.1.

After running this code and creating plots via "make .plots", you
should be able to view the plots in _plots/_PlotIndex.html.

