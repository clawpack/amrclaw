
.. _node-desc:

================
Node descriptors
================


.. _int-node:

------------------------
Integer node descriptors
------------------------

``node(levelptr,mptr)``
	Type: Pointer to grid
	
	Modified in: *birect*, *domain*, *grdfit*

	Description: Provides a pointer to the next grid on the same level, so that
	each level may be traversed as a linked list.  Value is set to null (global
	variable, currently =0) if mptr points to the last grid on the level.
	
``node(tempptr,mptr)``
	Type: Pointer to storage
	
	Target created in: *spest* (also used in plotclaw)
	
	Target destroyed in: *bufnst*

	Description: Provides the location of a temporary storage space where an error
	estimate is stored.  The allotted size is one word per grid cell per variable,
	including ghost cells.


``node(errptr,mptr)``
	Type: Pointer to storage
	
	Target created in: *errest*
	
	Target destroyed in: *errest*

	Description: Provides the location of a temporary storage space where another
	error estimate is stored, this estimate being on the grid *mptr* coarsened by
	a factor of 2.  The allotted size is one word per coarsened cell (including
	ghost cells) per variable (dependent and aux).

	Question: Does the coarsening correspond to the use of Richardson extrapolation?
	Yes.
	

``node(nestlevel,mptr)``
	Type: Integer
	
	Modified in: *domain*, *birect*, *grdfit*

	Description: Level containing the grid *mptr*.


``node(cfluxptr,mptr)``
	Type: Pointer to storage
	
	Target created in: *prepc*
	
	Target destroyed in: *putsp*

	Description: For each grid at level *coarseLevel* above the finest level,
	this provides the storage in a layer around each grid at level 
	*coarseLevel+1*.

	The space allotted is *5\*maxsp*, where *maxsp* is the number of level
	*coarseLevel* interfaces surrounding each of the finer grids.  This provides
	space for 5 fields to describe the interaction between grid *mptr* and the
	finer grid.  For later reference, *listsp(currentLevel)=maxsp*.
	
	Note that coarse fluxes themselves are not stored here.  They contribute
	instead to the space targeted by *ffluxptr* of each finer grid.
	
	To do: What are the five fields?  (See *setuse*: i, j, which side, kid, 
	location in kid's grid.)
	


``node(ffluxptr,mptr)``
	Type: Pointer to storage
	
	Target created in: *prepf*
	
	Target destroyed in: *putsp*

	Description: For each grid at level *fineLevel* below the coarsest level,
	this provides the storage location for fluxes in a layer around the grid,
	to be used in coarse-fine fixup.

	The space allotted is *2\*nvar\*lenbc+naux\*lenbc*, where *lenbc* is 2 times
	the number of boundary interfaces.  One space is for plus or minus fluxes, 
	and the other is for the coarse solution for wave fixing.
	
	
``node(store1,mptr)``
	Type: Pointer to storage
	
	Target created in: *ginit*, *gfixup*
	
	Target destroyed in: *gfixup*
	
	Description: Provides location in storage for the first copy of solution
	data.  The allotted size is (# of interior grid cells)*nvar.
	
	
``node(store2,mptr)``
	Type: Pointer to storage
	
	Target created in: *ginit*, *gfixup*
	
	Target destroyed in: *gfixup*

	Description: Provides location in storage for the second copy of solution
	data.  The allotted size is (# of interior grid cells)*nvar.


``node(ndilo,mptr)``
	Type: Integer
	
	Modified in: *domain*, *birect*, *grdfit*

	Description: Index of the leftmost interior cell in global index
	space.
	
	
``node(ndihi,mptr)``
	Type: Integer
	
	Modified in: *domain*, *birect*, *grdfit*

	Description: Index of the rightmost interior cell in global index
	space.


``node(ndjlo,mptr)``
	Type: Integer
	
	Modified in: *domain*, *birect*, *grdfit*

	Description: Index of the bottom-most interior cell in global index
	space.


``node(ndjhi,mptr)``
	Type: Integer
	
	Modified in: *domain*, *birect*, *grdfit*

	Description: Index of the top-most interior cell in global index
	space.


``node(storeaux,mptr)``
	Type: Pointer to storage
	
	Target created in: *ginit*, *gfixup*
	
	Target destroyed in: *gfixup*

	Description: Provides the location in storage designated for aux
	variables.  The allotted size is (# of interior cells)*naux.
	
	
``node(nextfree,mptr)``
    Type: Pointer to (free) node
    
    Description: Points to the next node on the free list.  Only
    relevant for free nodes, and currently uses the same space as
    ``node(tempptr,mptr)``.  See also: :ref:`nodget`, :ref:`putnod`.


.. _real-node:

---------------------
Real node descriptors
---------------------

``rnode(cornxlo,mptr)``
	Type: double precision
	
	Modified in: domain, grdfit, birect
	
	Description: Lower x-coordinate of the grid specified by mptr.
	

``rnode(cornylo,mptr)``
	Type: double precision
	
	Modified in: domain, grdfit, birect

	Description: Lower y-coordinate of the grid specified by mptr.


``rnode(cornxhi,mptr)``
	Type: double precision
	
	Modified in: domain, grdfit, birect

	Description: Upper x-coordinate of the grid specified by mptr.


``rnode(cornyhi,mptr)``
	Type: double precision
	
	Modified in: domain, grdfit, birect

	Description: Upper y-coordinate of the grid specified by mptr.		
	
	
``rnode(timemult,mptr)``
	Type: double precision
	
	Modified in: advanc, birect, ginit, grdfit, setgrd
	
	Description: Current time of the grid specified by mptr.
	
	Question: What is 'mult' in 'timemult'?
	"Multiple" of delta t
	But this actually stores the physical time