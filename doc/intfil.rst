
.. _intfil:

============
``intfil``
============

Signature: 
	``subroutine intfil(val, mi, mj, time, flaguse, nrowst, 
	ncolst, ilo, ihi, jlo, jhi, level, nvar, naux)``


Arguments: 
	``integer, intent(in)``:
		* ``mi``:   Horizontal dimension of the input region.
		* ``mj``:   Vertical dimension of the input region.
		* ``nrowst``: Starting ``i``-index of the patch,
		  in indices relative to ``val``.
		* ``ncolst``: Starting ``j``-index of the patch,
		  in indices relative to ``val``.
		* ``ilo``:  Lower horizontal index of the patch to fill,
		  in global indices.
		* ``ihi``:  Upper horizontal index of the patch to fill,
		  in global indices.
		* ``jlo``:  Lower vertical index of the patch to fill,
		  in global indices.
		* ``jhi``:  Upper vertical index of the patch to fill,
		  in global indices.
		* ``level``: Level of the input region.
		* ``nvar``: Number of solution variables.
		* ``naux``: Number of aux variables.
	``double precision, intent(inout)``:
		* ``val``:  Array of solution values on the patch being
		  filled.
	``double precision, intent(out)``:
		* ``flaguse(ilo:ihi,jlo:jhi)``: Indicates where the patch
		  is successfully filled.  ``flaguse`` is initialized to 
		  ``0.d0``, and set to ``1.d0`` in every location of the 
		  patch that gets filled.  It is also set to ``1.d0``
		  in all cells outside the computational domain, as these
		  will be filled in later using boundary data.
	``double precision, intent(in)``:
		* ``time``: Time on the patch being filled.
		
		
Description: 
	Attempts to fill a rectangular patch of ``val`` and ``aux``,
	which contain solution and aux values on an input region.  
	Values are needed at time ``time``, and at level ``level``.

	This routine only copies information from grids at level 
	``level`` into the patch; no spatial interpolation is done here.

	In  indices relative to ``val`` and ``aux``, the lower-left 
	corner of the patch is at ``(nrowst,ncolst)``.  In global 
	indices, its lower-left corner is at ``(ilo,jlo)``, and its 
	upper-right corner is at ``(ihi,jhi)``.

	This is a subprocess of ``filrecur``.