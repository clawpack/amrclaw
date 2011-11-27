
.. _setaux:

============
``setaux``
============

Signature: 
	``subroutine setaux(maxmx, maxmy, mbc, mx, my, xlower, ylower,
	dx, dy, maux, aux)``


Arguments: 
	``integer, intent(in):``
		* ``maxmx``: Maximum number of (interior) horizontal cells.
		* ``maxmy``: Maximum number of (interior) vertical cells.
		* ``mbc``:   Number of ghost cells in each direction.
		* ``mx``:    Actual number of horizontal cells.
		* ``my``:    Actual number of vertical cells.
		* ``maux``:  Number of aux variables.
	``double precision, intent(in):``
		* ``xlower``: Left bound of the input region.
		* ``ylower``: Bottom bound of the input region.
	    * ``dx``:	 Horizontal grid spacing.
		* ``dy``:	 Vertical grid spacing.
	``double precision, intent(out):``
		* ``aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)``:   Array 
		  of aux values, to be set by this routine.


Description: 
	User-supplied routine that initializes values of auxiliary
	variables.
	
**Question**:
	Why is this called at every time step?