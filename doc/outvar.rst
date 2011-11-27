
.. _outvar:

===========
``outvar``
===========

Signature:
    ``subroutine outvar(rect, mitot, mjtot, nvar, mptr, ng)``


Arguments:
    ``integer, intent(in)``:
        * ``mitot``: Horizontal dimension of grid, including ghost cells.
        * ``mjtot``: Vertical dimension of grid, including ghost cells.
        * ``nvar``: Number of solution variables.
        * ``mptr``: Pointer to grid being output.
    ``double precision, intent(in)``:
        * ``rect(mitot,mjtot,nvar)``: Solution values on the grid being
          output.


Description:
    Outputs the solution on a single grid, for later use by graphics routines.