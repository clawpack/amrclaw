
.. _bound:

============
``bound``
============

Signature:
    ``bound(time, nvar, ng, valbig, mitot, mjtot, mptr, aux, naux)``
	
Arguments:
    ``integer, intent(in)``:
        * ``nvar``: Number of solution variables.
        * ``ng``: Number of ghost cells in each direction (width of boundary region).
        * ``mitot``: Horizontal dimension of the grid, including ghost cells.
        * ``mjtot``: Vertical dimension of the grid, including ghost cells.
        * ``mptr``: Pointer to the descriptor of the input grid.
        * ``naux``: Number of aux variables.
    ``double precision, intent(in)``:
	    * ``time``: Time at which to fill the ghost cells.
	    * ``aux(mitot,mjtot,naux)``: Values of aux variables on the input grid.
    ``double precision, intent(inout)``:
        * ``valbig(mitot,mjtot,nvar)``: Solution values on the input grid.


Description:
    This routine sets the boundary values (ghost cells) for a given grid,
    specified by ``mptr``, at level ``level``.  It fills the values for
    a region ``ng`` cells wide all the way around the border, in 4
    rectangular strips.