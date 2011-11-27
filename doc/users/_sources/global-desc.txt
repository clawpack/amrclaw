
.. _global-desc:

=====================
Global descriptors
=====================

--------------
Parameters
--------------

All the below are of type ``integer, parameter``, and hard-coded
in ``call.i``.

``maxlv``
	Description: Maximum number of levels.
	

``maxgr``
	Description: Maximum number of grids, across all levels.
	

``maxcl``
	Description: Maximum number of clusters(?) used in regridding.
	

``max1d``
	Description: Maximum size of a grid along a single dimension.
	
	
``maxvar``
	Description: Maximum number of independent variables.
	
	
``maxaux``
	Description: Maximum number of aux variables.
	
	
``maxout``
	Description: Maximum number of times for output.
	
	
-----------------
Variables
-----------------

``ndfree``
    Description: Points to the head of the free list of nodes (grids).
    
    
``lfine``
    Description: Finest level currently in use.