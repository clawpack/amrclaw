
.. _putnod:

============
``putnod``
============


Signature:
    ``subroutine putnod(mptr)``
    

Arguments:
    ``integer, intent(in)``:
        * ``mptr``: Pointer to the grid being returned to the free list.
        

Description:
    Places the node ``mptr`` at the head of the free list.  This entire
    routine is simply::
        
        node(nextfree, mptr) = ndfree
        ndfree               = mptr