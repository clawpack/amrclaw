
.. _prepc:

=========
``prepc``
=========


Signature:
    ``subroutine prepc(level,nvar)``
    
    
Arguments:
    ``integer, intent(in)``:
        * ``level``: Coarse level, for which coarse flux
          descriptors are being prepared.
        * ``nvar``: Number of solution variables.
        
        
Description:
    After fine level ``level+1`` has been regridded, this routine
    allocates and fills (via :ref:`setuse`) the coarse flux 
    descriptors in ``node(cfluxptr,mcoarse)`` for each coarse grid
    ``mcoarse`` at level ``level``.