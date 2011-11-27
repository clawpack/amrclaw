
.. _outtre:

============
``outtre``
============

Signature: 
	``subroutine outtre(mlev, outgrd, nvar, naux)``


Arguments: 
	``integer, intent(in)``:
		* ``mlev``: Points to a grid on the coarsest level to 
		  output.  I.e., output begins with the level
		  ``node(nestlevel,mptr)``.
		* ``nvar``: Number of solution variables.
		* ``naux``: Number of aux variables.
	``logical, intent(in)``:
		* ``outgrd``: If ``outgrd=.true.``, then solution and
		  aux variables are output as well as grid descriptors.


Description:
	Output data from a subtree of the grid hierarchy.
