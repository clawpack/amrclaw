
.. _adj-struct:

=========================
Adjacency data structures
=========================


*node(adjxlo,mptr)*
	Type: Pointer to storage
	
	Target created in: *prepadj* (new routine, modeled after *prepf*)

	Target destroyed in: *putspadj* (new routine, modeled after *putsp*)
	
	Description: Target is an array of integers along the left boundary,
	corresponding to ghost cells.  Entries describe the relation of ghost
	cells to other grids at the same level, with possible values:

	* *madj*, pointer to grid at the same level the contains this ghost 
	  cell, if such a grid exists.
	
	* 0, if the ghost cell is interior, but does not overlap a grid of 
	  the same level.
	
	* -1, if the ghost cell is exterior to the domain.


*node(adjxhi,mptr)*
	Type: Pointer to storage

	Description: Same function as *node(adjxlo,mptr)*, but for ghost cells
	along the right boundary.
	

*node(adjylo,mptr)*
	Type: Pointer to integer data
	
	Description: Same function as *node(adjxlo,mptr)*, but for ghost cells
	along the bottom boundary.
	

*node(adjyhi,mptr)*
	Type: Pointer to integer data

	Description: Same function as *node(adjxlo,mptr)*, but for ghost cells
	along the top boundary.