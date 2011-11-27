
.. _prepadj:

=============
``prepadj``
=============

Signature:
    ``subroutine prepadj(level)``
    

Arguments:
    ``integer, intent(in)``:
        * ``level``: Level on which adjacency structures are
          being added.


Description: 
    Creates the adjacency data structures for all grids
    on level ``level``, and fills them by calling :ref:`setadj`.




