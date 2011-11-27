
.. _nodget:

===========
``nodget``
===========


Signature:
    ``integer function nodget(dummy)``


Arguments:
    ``double precision, intent(in)``:
        * ``dummy``: Unused dummy argument, provided simply to allow
          the function call.
          
          
Description:
    Returns and removes the head of the free list, stored in ``ndfree``.
    Updates ``ndfree`` accordingly.