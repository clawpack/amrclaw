
.. _armclaw_examples_acoustics_1d_adjoint:

Acoustics 1D -- Heterogeneous example (with adjoint flagging)
=============================================================

1D acoustics in a piecewise constant medium to illustrate reflection and
transmission at an interface.

The density and bulk modulus of the medium are specified in setrun.py,
along with a parameter determining what initial conditions to use
(see `qinit.f`).


Adjoint flagging
----------------

The adjoint method is used to flag cells needing refinement, as described in
the paper:

- Analysis and Performance Evaluation of Adjoint-Guided Adaptive Mesh
  Refinement for Linear Hyperbolic PDEs Using Clawpack, by
  B. N. Davis and R. J. LeVeque, 2018.
  `[link] <http://faculty.washington.edu/rjl/pubs/adjoint2018>`_

This example is similar to the problem in Example 3 of the paper.


Folder Organization
--------------------

- **adjoint:**

  Contains code to solve the adjoint problem.

  The output times specified in this directory should agree with those for the
  forward code.

Running the Code
--------------------

Go to the folder `adjoint` and run in a terminal::

    make new
    make .plots

The code will produce two new folders: _output and _plots. 
The first one contains all the output files, while the latter one contains
the plots and interactive visualization apps.

Go to the main folder `acoustics_2d_adjoint` and run in the terminal:

    make new
    make .plots

The code will produce two new folders: _output and _plots. 
The first one contains all the output files, while the latter one contains
the plots and interactive visualization apps.


Alternatively, to run first the adjoint and then the forward problem a
single script can be invoked.  
Run in the terminal::

    python run_adjoint_flagging.py

Running Variations
--------------------

In `setrun.py`, the following flags are currently set::

    flag_forward_adjoint = True
    flag_richardson_adjoint = False

This indicates that flagging will be done based only on the magnitude of the
inner product of the adjoint with the forward solution, using the tolerance
specified by the parameter `adjoint_flag_tolerance`.

Alternatively, the magnitude of the inner product of the adjoint with the
*estimated 1-step error* in the forward solution can be used to flag cells
for refinement, by setting::

    flag_forward_adjoint = False
    flag_richardson_adjoint = True

The same parameter `adjoint_flag_tolerance` is used to specify the level for
flagging.  This latter option is still under development, but is intended to
allow specifying a tolerance based on the level of global error hoped for in
the final solution.   See the paper cited above for more details.


This example can also be run with either the standard AMRClaw flagging based
on undivided differences (done in the `flag2refine` routine, where the
difference exceeds the tolerance `amrdata.flag2refine_tol` in absolute value)
or with Richardson extrapolation error flagging (in which case points are
flagged where the 1-step error is estimated to be above the tolerance
`amrdata.flag_richardson_tol`).

To revert to one of these older methods, in setrun.py` set::

    flag_forward_adjoint = False
    flag_richardson_adjoint = False

change one of these to `True`::

    amrdata.flag_richardson = False
    amrdata.flag2refine = False

and set the appropriate tolerance.
