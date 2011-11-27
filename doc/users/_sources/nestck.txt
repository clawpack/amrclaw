
.. _nestck:

============
``nestck``
============

Signature: 
	``logical function nestck(mnew, lbase, badpts, npts, numptc,
	icl, nclust, domflags, isize, jsize, nvar, naux)``


Arguments:
	``integer, intent(in)``:
		* ``mnew``:  Pointer to the grid being tested for proper containment.
		* ``lbase``: The finest level remaining fixed in this regridding.
		* ``nvar``:  Number of dependent variables.
		* ``naux``:  Number of aux variables.
	**Incomplete**:
		* ``badpts``
		* ``npts``
		* ``numptc``
		* ``icl``
		* ``nclust``
		* ``domflags``
		* ``isize``
		* ``jsize``

Description:
	Checks whether the potential grid ``mnew`` is completely contained
	in the (coarser) finest grid which stays fixed, at level ``lbase``.
	The projection algorithm will guarantee containment in all finer grids 
	between them.
