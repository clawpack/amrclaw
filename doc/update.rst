
.. _update:

============
``update``
============

Signature:
    ``subroutine update(level, nvar)``
    

Arguments:
    ``integer, intent(in)``:
        * ``level``: Level being updated.
        * ``nvar``: Number of solution variables.


Description:
    Updates coarse grids on level ``level`` based on fine grids
    at level ``level+1``.