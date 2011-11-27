
.. _cstore:

============
``cstore``
============

Signature: 
	``subroutine cstore(qc, nrow, ncol, nvar, qc1d, lenbc, naux, auxc, auxc1d)``


Arguments: 
	``integer, intent(in)``:
		* ``nrow``: Horizontal dimension of the coarse patch.
		* ``ncol``: Vertical dimension of the coarse patch.
		* ``nvar``: Number of solution variables.
		* ``lenbc``: Perimeter of the fine grid, in number of cells.
		* ``naux``: Number of aux variables.
	``double precision, intent(in)``:
		* ``qc(nrow,ncol,nvar)``: Solution values on coarse patch.
		* ``auxc(nrow,ncol,nvar)``: Aux variable values on coarse patch.
	``double precision, intent(out)``:
		* ``qc1d(lenbc)``: Stores coarse grid solution around perimeter
		  of fine grid.
		* ``auxc1d(lenbc)``: Stores coarse grid aux variables around perimeter
		  of fine grid.


Description: 
	Takes data (``qc`` and ``auxc``) from a coarse patch that perfectly surrounds
	a fine grid, and extracts the coarse data around the perimeter of the fine 
	grid.  These are stored in 1-dimensional arrays ``qc1d`` and ``auxc1d``.
	
	The perimeter of the fine grid is traversed in the following order::

		.          2
		.     __________
		.    |          |
		.  1 |          | 3		
		.    |__________|
		.          4

	save first interior cell of enlarged grid corresponding to
	fine grid bordering cell. note that since fine grid is smaller,
	the cell is one in. coarse (temporary) grid has no ghost cells
