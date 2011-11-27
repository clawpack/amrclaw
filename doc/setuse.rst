
.. _setuse:

==========
``setuse``
==========


Signature:
    ``subroutine setuse(listbc, maxsp, ispot, mkid, ilo, ihi,
    jlo, jhi, iclo, ichi, jclo, jchi, kflag)``


Arguments:
    ``integer, intent(in)``:
        * ``maxsp``:
        * ``mkid``:
        * ``ilo``:
        * ``ihi``:
        * ``jlo``:
        * ``jhi``:
        * ``iclo``:
        * ``ichi``:
        * ``jclo``:
        * ``jchi``:
        * ``kflag``:

    ``integer, intent(inout)``:
        * ``listbc(5,maxsp)``:
        * ``ispot``:



Description:
    set up boundary list for coarse grid, to be used by fluxsv. 
    loop around boundary of fine grids to do this.  each entry has
    i, j, side #, fine grid #, loc in fine grid list for fluxes.
    for example, side 1 of fine grid fixes side 3 of coarse grid,
    so coarse grid list will store the # 3.
    wrt coarse grid, the sides are::

        .          2
        .       1     3       that is, right edge of a coarse cell = 3
        .          4                    top  edge of a coarse cell = 2

    lkid is the index into the fine grid's saved fluxes.
    the fine grid will save all its fluxes all around its
    perimeter. lkid tells where the coarse grid should
    taking them from. (no ghost cells in this index, but 
    it is 1-based for indexing array, not - based for
    integer index of grid location).

    changed 11/11/08: spheredom for periodically mapped spherical
    grids. could affect top and bottom if fine grid touches
    edge of domain in y direction. if calling with spheredom
    (and not yperdom) then grid is NOT periodically mapped.
    need kflag to indicate spherically mapped now - otherwise
    cant tell the difference, dont skip appropropriate loops
    