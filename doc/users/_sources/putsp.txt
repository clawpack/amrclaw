
.. _putsp:

==========
``putsp``
==========

Signature: 
	``putsp(lbase,level,nvar,naux)``

Arguments:
	``integer, intent(in)``: 
		* ``lbase``: Base level of regridding; finer levels are being
		  regridded, but ``lbase`` is not.
		* ``level``: Level on which space is being reclaimed.
		* ``nvar``:  Number of independent variables.
		* ``naux``:  Number of dependent variables. 


Description: 
	Reclaims flux storage space in the main storage
	array for grids at level ``level``.  If ``level=lbase``, only
	the space at ``cfluxptr`` is reclaimed from each grid.  If 
	``level>lbase``, then the space at ``ffluxptr`` is reclaimed as well.