
.. _fixcapaq:

============
``fixcapaq``
============

Signature: 
	``subroutine fixcapaq(val, aux, mitot, mjtot, valc, auxc, mic, mjc,
	nvar, naux, levc)``


Arguments: 
    ``integer, intent(in):``
        * ``mitot``: Horizontal dimension of the fine grid.
        * ``mjtot``: Vertical dimension of the fine grid.
        * ``mic``: Horizontal dimension of the coarse patch.
        * ``mjc``: Vertical dimension of the coarse patch.
        * ``nvar``: Number of solution variables.
        * ``naux``: Number of aux variables.
        * ``levc``: Level of the coarse patch
    ``double precision, intent(in)``:
        * ``valc(mic,mjc,nvar)``: Solution values on the coarse patch.
        * ``auxc(mic,mjc,naux)``: Aux variable values on the coarse patch.
        * ``aux(mitot,mjtot,naux)``: Aux variable values on the fine grid.
    ``double precision, intent(inout)``:
        * ``val(mitot,mjtot,nvar)``: Solution values on the fine grid.
				

Description:
    After filling a new fine grid solution via linear interpolation,
    *kappa\*q* may need to be conserved rather than *q* in the presence
    of a capacity function.  This routine calculates the discrepancy in 
    *kappa\*q* and modifies *q* to account for it.

    The inputs are solution and aux data (``val``, ``aux``) from a fine 
    grid that has just been filled via linear interpolation, as well as 
    from a coarse patch (``valc``, ``auxc``) that covers the fine grid.