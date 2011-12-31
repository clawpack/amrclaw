""" 
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

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
    
    from clawutil.src.python import clawdata 
    
    assert claw_pkg.lower() == 'amrclaw',  "Expected claw_pkg = 'amrclaw'"

    ndim = 2
    rundata = clawdata.ClawRunData(claw_pkg, ndim)

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
    clawdata.ndim = ndim
    
    # Lower and upper edge of computational domain:
    clawdata.xlower = 0.0
    clawdata.xupper = 1.0
    
    clawdata.ylower = 0.0
    clawdata.yupper = 1.0
        

    # Number of grid cells:
    clawdata.mx = 50
    
    clawdata.my = 50
        

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 1

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 0
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.mcapa = 0
    
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        clawdata.output_ntimes = 10
        clawdata.output_tfinal = 2.0

    elif clawdata.output_style == 2:
        # Specify a list of output times.  
        clawdata.output_times =  [0.5, 1.0]   

    elif clawdata.output_style == 3:
        # Output every step_interval timesteps with a total of nsteps time steps:
        clawdata.output_step_interval = 1
        clawdata.output_nsteps = 5
        
    elif clawdata.output_style == 4:
        # Specify time interval for output up to tfinal:
        clawdata.output_time_interval = 0.2
        clawdata.output_tfinal = 2.0
    
    clawdata.output_format == 'ascii'      # 'ascii' or 'netcdf'
    clawdata.output_q_components = 'all'   # list of components to output
    clawdata.output_aux_components = 'all'
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

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = 1
    
    # Initial time step for variable dt.  
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.016
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used, and max to allow without 
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.0
    
    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000

    
    

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Transverse order for 2d or 3d (not used in 1d):
    clawdata.order_trans = 2
    
    # Number of waves in the Riemann solution:
    clawdata.mwaves = 1
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == mwaves
    clawdata.limiter = [3]
    
    # Source terms splitting:
    #   src_split == 0  => no source term (src routine never called)
    #   src_split == 1  => Godunov (1st order) splitting used, 
    #   src_split == 2  => Strang (2nd order) splitting used,  not recommended.
    clawdata.src_split = 0
    
    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.mbc = 2
    
    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity
    
    clawdata.bc_xlower = 'periodic'
    clawdata.bc_xupper = 'periodic'
                         
    clawdata.bc_ylower = 'periodic'
    clawdata.bc_yupper = 'periodic'
    

    # ---------------
    # AMR parameters:
    # ---------------


    # max number of refinement levels:
    amrlevel_max = 3

    # List of refinement ratios at each level (length at least mxnest+1)
    clawdata.refinement_ratio_x = [2,2,2]
    clawdata.refinement_ratio_y = [2,2,2]
    clawdata.refinement_ratio_t = [2,2,2]


    # Specify type of each aux variable in clawdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    clawdata.auxtype = []


    clawdata.flag_richardson = False  # don't use Richardson estimator
    clawdata.flag_gradient = True     # flag based on gradient of solution
    clawdata.flag_gradient_tol = 0.05  # used in default flag2refine subroutine
    clawdata.regrid_step_interval = 2  # how often to regrid (every kcheck steps)
    clawdata.regrid_buffer_width  = 3  # width of buffer zone around flagged points
    clawdata.verbosity_regrid = 3      # print regrid info up to this level

    # More AMR parameters can be set -- see the defaults in clawutil/??/clawdata.py

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
    
