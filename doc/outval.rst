
.. _outval:

============
``outval``
============

Signature: 
	``subroutine outval(val, nvar, mitot, mjtot, mptr, outgrd,
	naux,aux)``


Arguments: 
	``integer, intent(in)``:
		* ``nvar``:  Number of dependent variables.
		* ``mitot``: Number of horizontal cells, including ghost cells.
		* ``mjtot``: Number of vertical cells, including ghost cells.
		* ``mptr``:  Pointer to output grid.
		* ``naux``:  Number of aux variables.
	``logical, intent(in)``:
		* ``outgrd``: Only output if this is ``.true.``
	``double precision, intent(in)``:
		* ``val(mitot,mjtot,nvar)``: Solution values to print.
		* ``aux(mitot,mjtot,naux)``: Aux variables to print.


Description: 
	Prints the solution and aux variables to output; only prints values
	at interior (non-ghost) cells.

