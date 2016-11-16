Two-dimensional acoustics with radially symmetric initial data
==============================================================

Acoustics with radial symmetric initial conditions.  The solution should 
remain radially symmetric.  

First, run the code in the ‘adjoint’ subdirectory, and then run the code in
this directory. The code is set up to flag grid cells for refinement based 
on the inner product of the forward solution with the adjoint solution (or
on the inner product of the error of the forward solution with the adjoint
solution, if the Richardson flagging is selected).

