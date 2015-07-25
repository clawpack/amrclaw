
See docstring in run_with_restart.py

To test 2d version...

    cd test_advection_2d_square
    make .exe
    python ../run_with_restart.py
    Ctrl-C
    python ../run_with_restart.py

Similary in 3d.

The file test_advection_2d_square/check_twofiles.f should be moved to
amrclaw/src/2d eventually.

The file run_with_restart.py should be moved to clawutil/src/python/clawutil?
