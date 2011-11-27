
.. _intcopy:

============
``intcopy``
============

Signature: 
	``subroutine icall(val, mitot, mjtot, nvar, ilo, ihi, jlo, jhi, level,
	iputst, jputst)``


Arguments: 
	``integer, intent(in):``
		* ``mitot``: Horizontal dimension of the input patch.
		* ``mjtot``: Vertical dimension of the input patch.
		* ``nvar``: Number of solution variables.
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
		* ``val(mitot,mjtot,nvar)``: Solution values on the input patch.
				

Description: 
    Fills *only* solution data on an input patch with data from 
    intersecting grids at the same level.  To fill both solution and aux
    data, use :ref:`icall`.