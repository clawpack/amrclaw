
.. _saveqc:

============
``saveqc``
============

Signature: 
	``subroutine saveqc(level, nvar, naux)``


Arguments: 
	``integer, intent(in)``:
		* ``level``: Level of the fine grids, around which coarse fluxes
   		  are being saved.
		* ``nvar``: Number of solution variables.
		* ``naux``: Number of aux variables.


Description: 
	Loops over each grid on (fine) level ``level``.  For each such grid
	``mkid``:
	
	    #. Makes a coarsened, enlarged patch that extends one cell past 
	       ``mkid`` on each side.
	    #. Fills the coarsened patch via :ref:`icall` or ``preicall``.
	    #. Calls :ref:`cstore` to store values at the perimeter of the coarsened 
	       patch in the target of ``node(ffluxptr,mkid)``.