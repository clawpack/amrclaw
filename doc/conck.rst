
.. _conck:

==========
``conck``
==========

Signature:
    ``subroutine conck(level, nvar, time)``
    
    
Arguments:
    ``integer, intent(in)``:
        * ``level``: Level of the conservation check.  Most likely,
          ``level=1``, and conservation is being checked on the entire
          domain.
        * ``nvar``: Number of solution variables.  Not actually used.
    ``double precision, intent(in)``:
        * ``time``: Time of the conservation check.
        
    
Description:
    Outputs mass information on level ``level`` to standard output, so
    conservation may be inspected.  Primarily a debugging tool.  Assumes
    that grids don't overlap.