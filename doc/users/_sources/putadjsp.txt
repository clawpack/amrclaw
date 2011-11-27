
.. _putadjsp:

=============
``putadjsp``
=============


Signature:
    ``subroutine putadjsp(level)``
    
    
Arguments:
    ``integer, intent(in)``:
        * ``level``: Level from which adjacency structures are
          being reclaimed.
          

Description:
    Reclaim space for adjacency data structures for all grids
    on level ``level``.