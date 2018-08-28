- System requirement:
    - pgi compiler is installed; pgfortran and pgc++ are in system PATH (such that they can be called from anywhere)
    - CUDA 9.2 
    - CUDA_PATH environment variable is properly set to CUDA installation directory


- Installation (code setup):
    1. Install Clawpack by following the instruction from: http://www.clawpack.org/developers.html
    2. Add another git remote for the following clawpack sub-modules: amrclaw, geoclaw, clawutil, riemann. Then update the repo with git fetch and check out a different branch as below.
       Specifically, do the following for each sub-module:
       amrclaw:
            1) cd $CLAW/amrclaw
            2) git remote add xinshengqin https://github.com/xinshengqin/amrclaw.git
            3) git fetch xinshengqin
            4) git checkout merge_cudaclaw
       clawutil:
            1) cd $CLAW/clawutil
            2) git remote add xinshengqin https://github.com/xinshengqin/clawutil.git
            3) git fetch xinshengqin
            4) git checkout gpu_amr
       riemann:
            1) cd $CLAW/riemann
            2) git remote add xinshengqin https://github.com/xinshengqin/riemann.git
            3) git fetch xinshengqin
            4) git checkout gpu_amr
       geoclaw:
            1) cd $CLAW/geoclaw
            2) git remote add xinshengqin https://github.com/xinshengqin/geoclaw.git
            3) git fetch xinshengqin
            4) git checkout gpu

- Run examples
    Three examples can be found at:
        $CLAW/amrclaw/examples/GPU/acoustics_2d_heterogeneous_radial
        $CLAW/amrclaw/examples/GPU/acoustics_2d_radial
        $CLAW/amrclaw/examples/GPU/shallow_water_no_topo


