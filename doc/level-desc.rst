
.. _level-desc:

=================
Level descriptors
=================


.. _int-level:

-------------------------
Integer level descriptors
-------------------------
	
Each of the below descriptors is of type ``integer, dimension(maxlv)``.	

	
``icheck``
	Modified in: *stst1*, *tick*

	Description: Counts the number of time steps taken on the current
	level.  This determines when the level should have its error estimated,
	and finer levels regridded.
	
	Question: Does this count steps taken using this level's time step,
	or time steps of the coarsest level?


``intratx``
	Modified in: *amr2ez* (read from *amr2ez.data*)

	Description: Horizontal refinement ratio, used to obtain the next-finest
	level.


``intraty``
	Modified in: *amr2ez* (read from *amr2ez.data*)

	Description: Vertical refinement ratio, used to obtain the next-finest
	level.


``iregend``
	Modified in: *domain*, *grdfit*

	Description: The largest *i*-index used by a grid at this level.
	

``iregst``
	Modified in: *domain*, *grdfit*

	Description: The smallest *i*-index used by a grid at this level.


``iregsz``
	Modified in: *domain*, *restrt*, *restrt_hdf*

	Description: The horizontal width of the region (computational domain),
	measured in cells at this level.


``jregend``
	Modified in: *domain*, *grdfit*

	Description: The largest *j*-index used by a grid at this level.


``jregst``
	Modified in: *domain*, *grdfit*

	Description: The smallest *j*-index used by a grid at this level.
	
	
``jregsz``
	Modified in: *domain*, *restrt*, *restrt_hdf*

	Description: The vertical width of the region (computational domain),
	measured in cells at this level.	


``kratio``
	Modified in: *amr2ez* (set equal to intratx, and never changed)
	
	Description: Time refinement ratio, which is used to determine the time
	step used on the next-finest level.


``listsp``
	Modified in: *prepc*
	
	Description: Records ``maxsp`` for each level, which indicates the space
	(``5*maxsp``) alloted for coarse flux storage (target of ``cfluxptr``) to each 
	grid on that level.


``lstart``
	Modified in: *domain*, *stst1*, *setgrd*, *gfixup*

	Description: Pointer to the first grid on the level, where "first" refers
	to its location in the ``node`` data structure.


``newstl``
	Modified in: *regrid*, *grdfit*

	Description: Same function as ``lstart``, but used to build a new start list.
	Eventually copied into ``lstart`` in *setgrd*.  ``lstart`` cannot be overwritten
	until it has been used to interpolate values on the new grids, hence the new 
	space is needed.




.. _real-level:

----------------------
Real level descriptors
----------------------

Each of the descriptors below is of type ``double precision, dimension(maxlv)``.


``hxposs``
	Modified in: *amr2ez*, *stst1*, *restrt*, *restrt_hdf*
	
	Description: Records *hx*, the horizontal cell width.


``hyposs``
	Modified in: *amr2ez*, *stst1*, *restrt*, *restrt_hdf*

	Description: Records *hy*, the vertical cell width.


``possk``
	Modified in: *amr2ez*, *stst1*, *restrt*, *restrt_hdf*

	Description: Records *k*, the length of the time step.