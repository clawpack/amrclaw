
.. _reclam:

============
``reclam``
============

Signature: 
	``subroutine reclam(index, nwords)``


Arguments: 
	``integer, intent(in):``
		* ``index``:  Starting location of space to be freed
		  in ``alloc``.
		* ``nwords``: Length of space to be freed.


Description: 
	Return ``nwords`` of space, beginning at location ``index``,
	to the free list.

