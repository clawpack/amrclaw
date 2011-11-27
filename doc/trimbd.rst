
.. _trimbd:

============
``trimbd``
============

Signature: 
	``subroutine trimbd(used, nrow, ncol, set, il, ir, jb, jt)``


Arguments: 
	``integer, intent(in)``:
		* ``nrow``: Horizontal dimension of ``used``.
		* ``ncol``: Vertical dimension of ``used``.
	``double precision, intent(in)``:
		* ``used(nrow,ncol)``: Array of flags indicating whether 
   		  solution values have been filled in a patch using ``filrecur``.
		  Value is ``0.d0`` at unset points, and ``1.d0`` at set points.
	``integer, intent(out)``:
		* ``il``: If ``used`` is not completely set, indicates the
		  lower horizontal index of the smallest rectangle containing
		  the unset points.
		* ``ir``: Like ``il``, but indicates the upper horizontal index.
		* ``jb``: Like ``il``, but indicates the lower vertical index.
		* ``jt``: Like ``il``, but indicates the upper vertical index.
	``logical, intent(out)``:
		* ``set``: Returns ``.true.`` if all elements of the ``used``
		  array are set.  Returns ``.false.`` otherwise.


Description: 
	If the ``used`` array is completely set (``=1.d0``) then this 
	routine returns	``set=.true.``.  Otherwise it returns 
	``set=.false.``, and specifies the smallest rectangle containing
	all unset points, which has lower-left corner ``(il,jb)`` and
	upper-right corner ``(ir,jt)``.