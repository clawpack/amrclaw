# setrun file for 2D acoustics in a homogeneous domain

def setrun(claw_pkg='classic'):
    from clawpack.clawutil import data

    # 2D general data object
    rundata = data.ClawRunData(claw_pkg, 1)    # 2 = number of dimensions

    # Problem-specific data
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('dim',2,  'ndim = space dimensions (2 = cylindrical symmetry, 3 = spherical)')
    probdata.add_param('rho',1.,  'density of medium')
    probdata.add_param('bulk',4.,  'bulk modulus')


    # General Clawpack control data
    clawdata = rundata.clawdata    # Initialized when rundata is created

    # Size of computational domain
    clawdata.lower = 0.0
    clawdata.upper = 1.5

    # Grid dimensions
    clawdata.num_cells = 2000;

    # Size of system
    clawdata.num_eqn = 2
    clawdata.num_aux = 0    # Original 4.6 code had maux = 2.  Typo?  Nothing used aux.
    clawdata.capa_index = 0    # Obviously

    # Output control
    clawdata.output_style = 1
    clawdata.num_output_times = 30
    clawdata.tfinal = 0.27

    # Time stepping
    clawdata.dt_initial = 0.5
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.
    clawdata.steps_max = 5000
    clawdata.dt_variable = True

    # Details of the numerical method
    clawdata.order = 2
    clawdata.verbosity = 1
    clawdata.source_split = 1
    clawdata.dimensional_split = 0
    clawdata.use_fwaves = False

    # Waves and limiter
    clawdata.num_waves = 2
    clawdata.limiter = [2, 2]

    # Boundary conditions
    # Zero-order extrapolation everywhere
    clawdata.bc_lower = [3]
    clawdata.bc_upper = [1]

    return rundata


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
	rundata = setrun()

    rundata.write()
