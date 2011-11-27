
.. _icall:

============
``icall``
============

Signature: 
	``subroutine icall(val, aux, nrow, ncol, nvar, naux, ilo, ihi, jlo,
	jhi, level, iputst, jputst)``


Arguments: 
	``integer, intent(in):``
		* ``nrow``: Horizontal dimension of the input patch.
		* ``ncol``: Vertical dimension of the input patch.
		* ``nvar``: Number of solution variables.
		* ``naux``: Number of aux variables.
		* ``ilo``: Lower horizontal index of the patch, in global indices.
		* ``ihi``: Upper horizontal index of the patch, in global indices.
		* ``jlo``: Lower vertical index of the patch, in global indices.
		* ``jhi``: Upper vertical index of the patch, in global indices.
		* ``level``: Level of the patch
		* ``iputst``: Lower horizontal index of the patch region to fill,
		  in indices local to the patch.
		* ``jputst``: Lower vertical index of the patch region to fill,
		  in indices local to the patch.
	``double precision, intent(inout)``:
		* ``val(nrow,ncol,nvar)``: Solution values on the input patch.
		* ``aux(nrow,ncol,nvar)``: Values of aux variables on the input patch.
				

Description: 
	Fills solution and aux data on an input patch with data from 
	intersecting grids at the same level.


**Notes:**
	Currently ``sticksout`` is passed as an extra argument when this
	routine is called from :ref:`saveqc`.