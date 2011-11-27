
.. _check:

===========
``check``
===========

Signature:
    ``subroutine check(nsteps, time, nvar, naux)``
    

Arguments:
    ``integer, intent(in)``:
        * ``nsteps``: Number of steps taken on the coarse grid since
          the start.
        * ``nvar``: Number of solution variables.
        * ``naux``: Number of aux variables.
    ``double precision, intent(in)``:
        * ``time``: Time of the checkpoint.


Description:
    Creates a checkpoint, from which the simulation may be restarted.
    Can only be called at the end of coarse grid cycles.