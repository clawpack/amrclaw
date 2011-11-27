
.. _basic:

==========
``basic``
==========

Signature:
    ``subroutine basic(time, lst, end)``
    

Arguments:
    ``integer, intent(in)``:
        * ``lst``: Coarsest level for output of tree structure.
        * ``lend``: Finest level for output of tree structure.
    ``double precision, intent(in)``:
        * ``time``: Time of output.
        
        
Description:
    Outputs basic information needed by the other graphics
    output routines (:ref:`valout`) at the given ``time``.
    Writes the entire level list, from level ``1`` to ``lfine``,
    and the tree structure from level ``lst`` to ``lend``.
