""" 
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

import os
from numpy import linspace, arange  # often used to set input arrays
from clawutil import clawdata 

#------------------------------
def setrun(claw_pkg='amrclaw'):
#------------------------------
    
    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "amrclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """ 
    
    
    assert claw_pkg.lower() == 'amrclaw',  "Expected claw_pkg = 'amrclaw'"

    num_dim = 2
    rundata = clawdata.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('u',     0.5,  'ubar advection velocity')
    probdata.add_param('v',     1.0,  'vbar advection velocity')
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim
    
    # Lower and upper edge of computational domain:
    clawdata.lower = [0., 1.]   # [xlower, ylower]
    clawdata.upper = [0., 1.]   # [xupper, yupper]

    ## Or is it clearer to write: ??
    #clawdata.lower[0] = 0.  # xlower
    #clawdata.upper[0] = 1.  # xupper
    #
    #clawdata.lower[1] = 0.  # ylower
    #clawdata.upper[1] = 1.  # yupper
    

    # Number of grid cells:
    clawdata.num_cells = [50,50]  # [mx,my]
    

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 1

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 0
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 0
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0                        
    
    # If restarting, t0 above should be from original run
    clawdata.restart = False                # True to restart from prior results
    clawdata.restart_file = 'restart.data'  # File to use for restart data
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

#===============================================
# Get rid of this old style !!!  See below.
#
#   clawdata.output_style = 1
#
#   if clawdata.output_style==1:
#       # Output ntimes frames at equally spaced times up to tfinal:
#       clawdata.num_output_times = 4
#       clawdata.tfinal = 2.0
#
#   elif clawdata.output_style == 2:
#       # Specify a list of output times.  
#       clawdata.output_times =  [0.5, 1.0, 1.5, 2.0]   
#
#   elif clawdata.output_style == 3:
#       # Output every step_interval timesteps over total_steps timesteps:
#       clawdata.output_step_interval = 2
#       clawdata.total_steps = 6
#       
#   elif clawdata.output_style == 4:
#       # Specify time interval for output up to tfinal:
#       clawdata.output_time_interval = 0.5
#       clawdata.tfinal = 2.0
#===============================================

    #===============================================
    # Instead always allow both a list of output times and
    # a list of output steps...  
    # Output will be produced after each specified step on the 
    # coarse grid (if any specified) 
    # AND at each specified time (by shortening time steps as needed).

    clawdata.output_times = linspace(0., 2., 5)
                     # or = [0, .5, 1., 1.5, 2.]
                     # or = arange(0,2.1,0.5)
                     # or = []  # to not specify output times
                     # etc.

    clawdata.output_steps = range(0,7,2)
                     # or = [0, 2, 4, 6]
                     # or = []  # to not specify output steps
                     # etc.

    # Notes: 
    # All previous output styles can be specified by one of the examples above.
    # If output_times == [] and output_steps == [] then no steps will be taken.
    # Output at t=t0 was always done previously.  Should it be?
    # Normally either output_times == [] or output_steps == [] and the other
    # determines when output is done, but no need to require this.
    

    clawdata.output_format == 'ascii'      # 'ascii' or 'netcdf' or 'binary'

    clawdata.output_q_components = 'all'   # could be list such as [0,2]
    clawdata.output_aux_components = []
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0
    

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==False: fixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True
    
    # Initial time step for variable dt.  
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 0.016
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used 
    clawdata.cfl_desired = 0.9
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.0
    
    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 1000000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 0
    
    # For unsplit method, transverse_waves can be 
    #  0 ==> donor cell (only normal solver used)
    #  1 ==> corner transport of waves
    #  2 ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2
    
    
    # Number of waves in the Riemann solution:
    clawdata.num_waves = 1
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 = no limiter (Lax-Wendroff)
    #   1 = minmod
    #   2 = superbee
    #   3 = MC limiter
    #   4 = van Leer
    clawdata.limiter = [3]
    
    clawdata.fwave = False    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0  => no source term (src routine never called)
    #   src_split == 1  => Godunov (1st order) splitting used, 
    #   src_split == 2  => Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 0
    
    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2
    
    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity
    
    clawdata.bc_lower[0] = 'periodic'   # at xlower
    clawdata.bc_upper[0] = 'periodic'   # at xupper

    clawdata.bc_lower[1] = 'periodic'   # at ylower
    clawdata.bc_upper[1] = 'periodic'   # at yupper
                         
    

    # ---------------
    # AMR parameters:
    # ---------------


    # max number of refinement levels:
    clawdata.amr_levels_max = 3

    # List of refinement ratios at each level (length at least amr_level_max-1)
    clawdata.refinement_ratio_x = [2,2]
    clawdata.refinement_ratio_y = [2,2]
    clawdata.refinement_ratio_t = [2,2]
    clawdata.refinement_ratio_t = [2,2]

    # Instead of setting refinement ratios in t, these can be chosen
    # automatically if this is implemented:
    clawdata.variable_dt_refinement_ratios = True
    # Currently only available in GeoClaw, where it's needed for special
    # case of tsunami modeling.  Is it useful more generally??


    # Specify type of each aux variable in clawdata.auxtype.
    # This must be a list of length num_aux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    clawdata.aux_type = []


    # Flag for refinement based on Richardson error estimater:
    clawdata.flag_richardson = False    # use Richardson?
    clawdata.flag_richardson_tol = 0.1  # Richardson tolerance
    
    # Flag for refinement based routine flag2refine:
    clawdata.flag2refine = True      # use this?
    clawdata.flag2refine_tol = 0.05  # tolerance used in this routine
    # User can modify flag2refine to change the criterion for flagging.
    # Default: check maximum absolute difference of first component of q
    # between a cell and each of its neighbors.

    # steps to take on each level L between regriddings of level L+1:
    clawdata.regrid_interval = 2       

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    clawdata.regrid_buffer_width  = 2  

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    clawdata.clustering_cutoff = 0.7  

    # print info about each regridding up to this level:
    clawdata.verbosity_regrid = 3      


    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 2

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at the final time.
        pass

    elif clawdata.checkpt_style == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = [1., 2.]

    elif clawdata.checkpt_style == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 100

    elif clawdata.checkpt_style == 4:
        # Checkpoint every checkpt_time_interval time units
        clawdata.checkpt_time_interval = 1.


    #  ----- For developers ----- 
    # Toggle debugging print statements:
    clawdata.dprint = True       # print domain flags
    clawdata.eprint = True       # print err est flags
    clawdata.edebug = True       # even more err est flags
    clawdata.gprint = False      # grid bisection/clustering
    clawdata.nprint = False      # proper nesting output
    clawdata.pprint = False      # proj. of tagged points
    clawdata.rprint = False      # print regridding summary
    clawdata.sprint = False      # space/memory output
    clawdata.tprint = True       # time step reporting each level
    clawdata.uprint = False      # update/upbnd reporting
    
    return rundata

    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
        rundata = setrun()

    rundata.write()
    
