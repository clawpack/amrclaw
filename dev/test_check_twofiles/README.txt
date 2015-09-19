
See docstring in run_with_restart.py

To test 2d version...

    cd test_advection_2d_square
    make .exe
    python ../run_with_restart.py
    Ctrl-C
    python ../run_with_restart.py
    [repeat this process as many times as desired until it finishes]

Then check the following:
    run_output.txt  contains stdout output from each run and info about restarts
    _output/fort.gauge  has final set of gauge data from last execution
    _output/fort.gauge_DATETIME  files contain fort.gauge from earlier
            
All fort.gauge* files need to be catenated together, but there may be some
overlapping times (after last checkpoint and before code died).

Similary in 3d.

The file test_advection_2d_square/check_twofiles.f should be moved to
amrclaw/src/2d eventually.

The file run_with_restart.py should be moved to clawutil/src/python/clawutil?
