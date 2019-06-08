
.. _armclaw_examples_acoustics_2d_adjoint:

Acoustics 2D -- Heterogeneous example (with adjoint flagging)
=============================================================

2D acoustics in a piecewise constant medium to illustrate reflection and
transmission at an interface.

The density and bulk modulus of the medium are specified in setrun.py.

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

  The output times specified in this directory should to at least as
  far out in time as the forward solution is desired, with sufficiently
  dense outputs to properly capture the evolving adjoint solution.

Running the Code
--------------------

Go to the folder `adjoint` and run in a terminal::

    make new
    make .plots

The code will produce two new folders: _output and _plots. 
The first one contains all the output files, while the latter one contains the plots and interactive visualization apps.

Then return to this directory and 

    make new
    make .plots

to run the forward solution and make plots that also show the inner product
of the forward and adjoint solution on levels 1 and 2 (not on level 3 since 
there is no need to flag further in this 3-level run).

Alternatively, to run first the adjoint and then the forward problem a
single script can be invoked.  
Run in the terminal::

    python run_adjoint_flagging.py

Running Variations
--------------------

In `setrun.py`, the following flags are currently set (in various places)::

    adjointdata.use_adjoint = True

    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False
    amrdata.flag_richardson_tol = 0.1
    
    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True
    rundata.amrdata.flag2refine_tol = 0.04

    # time period of interest:
    adjointdata.t1 = rundata.clawdata.t0
    adjointdata.t2 = rundata.clawdata.tfinal

Setting `adjointdata.use_adjoint` to `False` will go back to using standard
flagging based on the magnitude of undivided differences or an estimate of
the one-step error.

With adjoint flagging, either `amrdata.flag2refine` or
`amrdata.flag_richardson` can be set to `True` and the associated tolerance
is then applied to the inner product of the forward solution with the adjoint
solution (with `flag2refine`) or with the adjoint 1-step error estimate
(with `flag_richardson`).  
This latter option is still under development, but is intended to
allow specifying a tolerance based on the level of global error hoped for in
the final solution.   See the paper cited above for more details.

The time period of interest can be changed to some subset of the full run
time.  Try changing this to see how the AMR adapts to only capture waves
reaching Gauge 0 over the specified time period.

Note that the location of interest is specified in `adjoint/qinit.f90`, where
the functional used as initial data (at the final time) in the adjoint
problem is set to be a small rectangle around this location.


