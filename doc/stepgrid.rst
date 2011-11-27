
.. _stepgrid:

=============
``stepgrid``
=============

Signature:
    ``subroutine stepgrid(q, fm, fp, gm, gp, mitot, mjtot, mbc, dt, 
    dtnew, dx, dy, nvar, xlow, ylow, time, mptr, maux, aux)``
    

Arguments:
    ``integer, intent(in)``:
        * ``mbc``: Number of ghost cells.
        * ``mitot``: Horizontal dimension of grid ``mptr``, including
          ghost cells.
        * ``mjtot``: Vertical dimension of grid ``mptr``, including
          ghost cells.
        * ``nvar``: Number of solution variables.
        * ``mptr``: Pointer to grid being stepped.
        * ``maux``: Number of aux variables.
    ``double precision, intent(in)``:
        * ``dt``: Incoming time step.
        * ``dx``: Horizontal cell width of grid ``mptr``.
        * ``dy``: Vertical cell width of grid ``mptr``.
        * ``xlow``: Lower *x*-index of grid, including ghost cells.
        * ``ylow``: Lower *y*-index of grid, including ghost cells.
    ``double precision, intent(inout)``:
        * ``q(mitot,mjtot,nvar)``: Solution values on the grid, to be
          overwritten.
        * ``aux(mitot,mjtot,maux)``: Aux variable values on the grid.
    ``double precision, intent(out)``:
        * ``fm(mitot,mjtot,nvar)``: Fluxes left of cell edges.
        * ``fp(mitot,mjtot,nvar)``: Fluxes right of cell edges.
        * ``gm(mitot,mjtot,nvar)``: Fluxes below cell edges.
        * ``gp(mitot,mjtot,nvar)``: Fluxes above cell edges
        * ``dtnew``: Suggested new time step for this grid's solution.


Description:
    Takes a time step on the grid ``mptr``.