
.. _bc2amr:

============
``bc2amr``
============

Signature: 
	``subroutine bc2amr(val, aux, nrow, ncol, meqn, naux, hx, hy,
	level, time, xleft, xright, ybot, ytop, xlower, ylower, xupper,
	yupper, xperiodic, yperiodic, spheredom)``


Arguments: 
	``integer, intent(in)``:
		* ``nrow``: Horizontal dimension of the input patch.
		* ``ncol``: Vertical dimension of the input patch.
		* ``meqn``: Number of solution variables.
		* ``naux``: Number of aux variables.
		* ``level``: Level of the input patch.
	``double precision, intent(in)``:
		* ``time``:   Time at which boundary conditions are being applied.
		* ``xleft``:  Left bound of the input patch.
		* ``xright``: Right bound of the input patch.
		* ``ybot``:   Bottom bound of the input patch.
		* ``ytop``:   Top bound of the input patch.
		* ``xlower``: Left bound of the physical domain.
		* ``xupper``: Right bound of the physical domain.
		* ``ylower``: Bottom bound of the physical domain.
		* ``yupper``: Top bound of the physical domain.
		* ``aux(nrow,ncol,meqn)``: Aux variables on the input patch.
	``logical, intent(in)``:
		* ``xperiodic``: Value is ``.true.`` if the physical domain
		  is periodic in *x*.
		* ``yperiodic``: Value is ``.true.`` if the physical domain
		  is periodic in *y*.
		* ``spheredom``: Value is ``.true.`` if the physical domain
		  is spherically symmetric.
	``double precision, intent(inout)``:
		* ``val(nrow,ncol,meqn)``: Solution values on the input patch.


Description: 
	Takes a grid patch with mesh widths ``hx`` and ``hy``, of 
	dimensions ``nrow`` by ``ncol``, and uses the boundary conditions
	to set the solution values on any piece	of the patch which extends 
	outside the physical domain.
	
	User-specified boundary conditions must be included in this routine.