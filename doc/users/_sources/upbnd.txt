
.. _upbnd:

==========
``upbnd``
==========

Signature:
    ``subroutine upbnd(listbc, val, nvar, mitot, mjtot, maxsp,
    iused, mptr)``
    
    
Arguments:
    ``integer, intent(in)``:
        * ``listbc(5,maxsp)``: Coarse boundary lists for the input grid.
        * ``nvar``: Number of solution variables.
        * ``mitot``: Horizontal dimension of the input grid, including
          ghost cells.
        * ``mjtot``: Vertical dimension of the input grid, including
          ghost cells.
        * ``maxsp``: Space needed for each field in the coarse boundary
          lists.
        * ``mptr``: Pointer to the coarse grid being corrected.
    ``integer, intent(out)``:
        * ``iused(mitot,mjtot)``: Used to indicate where flux updates
          have been performed.
    ``double precision, intent(inout)``:
        * ``val(mitot,mjtot,nvar)``: Array of solution values on grid
          ``mptr``.


Description:
    Corrects the coarse grid ``mptr`` with the flux differences stored
    with each of the fine grids. 