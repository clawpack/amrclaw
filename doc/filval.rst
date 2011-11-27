
.. _filval:

============
``filval``
============

Signature: 
	``subroutine filval(val, mitot, mjtot, hx, hy, lev, time, valc, auxc,
	mic, mjc, xleft, xright, ybot, ytop, nvar, mptr, ilo, ihi, jlo, jhi,
	aux, naux, locflip)``


Arguments:
    ``integer, intent(in)``:
        * ``mitot``: Horizontal dimension of the fine grid, including ghost cells.
        * ``mjtot``: Vertical dimension of the fine grid, including ghost cells.
        * ``lev``: Level of the fine grid.
        * ``nvar``: Number of solution variables.
        * ``mptr``: Pointer to the fine grid being filled.
        * ``naux``: Number of aux variables.
        * ``ilo``: Lower horizontal index of the fine grid, including ghost cells, in global indices.
        * ``ihi``: Upper horizontal index of the fine grid, including ghost cells, in global indices.
        * ``jlo``: Lower vertical index of the fine grid, including ghost cells, in global indices.
        * ``jhi``: Upper vertical index of the fine grid, including ghost cells, in global indices.
        * ``locflip``: Location of perimeter storage for fine grid, for use in the periodic case.
        
    ``double precision, intent(in)``:
        * ``hx``: Horizontal cell width on the fine grid.
        * ``hy``: Vertical cell width of the fine grid.
        * ``time``: Physical time.
        * ``xleft``: Left coordinate of the fine grid, ghost cells excluded.
        * ``xright``: Right coordinate of the fine grid, ghost cells excluded.
        * ``ybot``: Bottom coordinate of the fine grid, ghost cells excluded.
        * ``ytop``: Top coordinate of the fine grid, ghost cells excluded.
        
    ``double precision, intent(inout)``:
        * ``val(mitot,mjtot,nvar)``: Solution values on the fine grid, with
          space for ghost cells.
        * ``aux(nrow,ncol,nvar)``: Values of aux variables on the fine grid,
          with space for ghost cells.
        
    ``double precision, intent(out)``:
	    * ``valc(mic,mjc,nvar)``: Solution values on the coarse patch.
	    * ``auxc(mic,mjc,naux)``: Aux variable values on the coarse patch.
				

Description:
    Fills in solution values on the interior of a fine grid.  This is done 
    by filling a coarse patch that extends one coarse cell beyond the interior
    of the fine grid on each side.  Then solution values on the fine grid may
    be linearly interpolated.

    Afterwards, solution values on the fine grid are replaced by values copied
    from other fine grids, wherever possible.