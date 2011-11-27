
.. _filrecur:

============
``filrecur``
============

Signature: 
	``recursive subroutine filrecur(level, nvar, valbig, aux,
	naux, time, mitot, mjtot, nrowst, ncolst, ilo, ihi, jlo, 
	jhi)``


Arguments: 
	``integer, intent(in)``:
		* ``level``: Level of the input region.
		* ``nvar``: Number of solution variables.
		* ``naux``: Number of aux variables.
		* ``mitot``: Horizontal dimension of the input region.
		* ``mjtot``: Vertical dimension of the input region.
		* ``nrowst``: Starting ``i``-index of the patch,
		  in indices relative to ``valbig``.
		* ``ncolst``: Starting ``j``-index of the patch,
		  in indices relative to ``valbig``.
		* ``ilo``:  Lower horizontal index of the patch to fill,
		  in global indices.
		* ``ihi``:  Upper horizontal index of the patch to fill,
		  in global indices.
		* ``jlo``:  Lower vertical index of the patch to fill,
		  in global indices.
		* ``jhi``:  Upper vertical index of the patch to fill,
		  in global indices.
	``double precision, intent(in)``:
		* ``time``: Time at which values are needed.
	``double precision, intent(inout)``:
		* ``valbig(mitot,mjtot,nvar)``: Solution values on the
		  input region.
		* ``aux(mitot,mjtot,naux)``: Aux variable values on the
		  input region.



Description:
	Fills in a rectangular patch of ``valbig`` and ``aux``, which
	contain solution and aux values on an input region.  Values 
	are needed at time ``time``, and at level ``level``.

	In  indices relative to ``valbig`` and ``aux``, the lower-left 
	corner of the patch is at ``(nrowst,ncolst)``.  In global 
	indices, its lower-left corner is at ``(ilo,jlo)``, and its 
	upper-right corner is at ``(ihi,jhi)``.

	First, the patch is filled with values obtainable from the 
	level ``level``	grids.   If any values are left unfilled, 
	the remaining rectangle of unfilled values is enlarged by 1
	(for later linear interpolation), and remaining values are
	recursively obtained from coarser levels.