
.. _outmsh:

============
``outmsh``
============

Signature: 
	``subroutine outmsh(mptr,outgrd,nvar,naux)``


Arguments: 
	``integer, intent(in)``:
		* ``mptr``: Pointer to the grid descriptor being output.
		* ``nvar``: Number of dependent variables.
		* ``naux``: Number of aux variables.
	``logical, intent(in)``:
		* ``outgrd``: If ``outgrd=.true.``, then solution values
		  and aux variables on the grid are output as well.


Description: 
	Outputs the grid descriptor, and optionally the values on the 
	grid, for a single grid referenced by ``mptr``.  (See 
	``outtre`` for outputting a full subtree.)
