
.. _fluxad:

===========
``fluxad``
===========

Signature:
    ``subroutine fluxad(xfluxm, xfluxp, yfluxm, yfluxp, svdflx,
    mptr, mitot, mjtot, nvar, lenbc, lratiox, lratioy, ng, dtf,
    dx, dy)``
    

Arguments:
    ``integer, intent(in)``:
        * ``mptr``:
        * ``mitot``:
        * ``mjtot``:
        * ``nvar``:
        


Description:
    save fine grid fluxes  at the border of the grid, for fixing
    up the adjacent coarse cells. at each edge of the grid, only
    save the plus or minus fluxes, as necessary. For ex., on
    left edge of fine grid, it is the minus xfluxes that modify the
    coarse cell.


Incomplete!