

.. _index:

***************
AMRClaw
***************

The $CLAW/amrclaw directory contains a version of Clawpack that uses
adaptive mesh refinement (AMR).

See :ref:`setrun_amrclaw` for a description of the run-time parameters
needed by AMRClaw and how to set them.

The refinement strategy and criteria are described in :ref:`amr_strategy`.

Several examples can be found from :ref:`apps`, look for ones with an *amr*
subdirectory.

The sections below detail many components of the library $CLAW/amrclaw/2d.



-----------------------------------
Global parameters and descriptors
-----------------------------------

This describes most of the entries in the common block, call.i.

.. toctree::
   :maxdepth: 2

   global-desc
   level-desc
   node-desc

-------------
Main program
-------------

.. toctree::
   :maxdepth: 2
   
   amr2ez


----------
Functions
----------

.. toctree::
   :maxdepth: 2

   igetsp
   nestck
   nodget



---------------------
Selected subroutines
---------------------

.. toctree::
   :maxdepth: 2

   basic
   bc2amr
   bound
   check
   cleanup
   conck
   cstore
   icall
   intcopy
   intfil
   filrecur
   filval
   fixcapaq
   flglvl
   outmsh
   outtre
   outval
   outvar
   prepc
   putnod
   putsp
   reclam
   saveqc
   setaux
   setuse
   stepgrid
   trimbd
   upbnd
   update
   valout
