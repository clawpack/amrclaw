
.. _flglvl:

============
``flglvl``
============

Signature: 
	``subroutine flglvl(nvar, naux, lcheck, nxypts, index, lbase,
	ldom2, npts, t0)``


Arguments: 
	``integer, intent(in)``:
		* ``nvar``:   Number of dependent variables.
		* ``naux``:   Number of aux variables.
		* ``lcheck``: Level to be flagged.
		* ``lbase``:  Finest level remaining fixed during regridding.
	``integer, intent(out)``:
		* ``nxypts``: Total number of flagged points.
		* ``index``:  Starting index in ``alloc`` of the flagged points
		  (which occupy ``2\*nxypts`` locations).
	**Incomplete**:
		* ``ldom2``
		* ``t0``
		

Description:
	Controls error estimation and/or flagging of bad points for the
	input level ``lcheck``.