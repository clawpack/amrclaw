
First pass at a restart example...

For full run and to create checkpoint files:
    make .plots 

Create new output directory for restart and copy over checkpoint files:
    mkdir _output_restart
    cp _output/fort*ch* _output_restart/

To test restart:
    make .plots -f Makefile_restart  

To compare plots:
    $CLAW/clawutil/src/python/clawutil/imagediff.py _plots _plots_restart

To make side-by-side plots (after running both codes, use plots not .plots):
    make plots -f Makefile_compare
