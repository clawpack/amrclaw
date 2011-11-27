
.. _cleanup:

============
``cleanup``
============

Signature:
    ``subroutine cleanup(nvar, naux)``
    

Arguments:
    ``integer, intent(in)``:
        * ``nvar``: Number of solution variables.
        * ``naux``: Number of aux variables.
        

Description:
    Final subroutine called by :ref:`amr2ez`.  Reclaims all remaining
    storage space from grids, allowing ``amr2ez`` to raise a warning
    if memory has not been properly released.
