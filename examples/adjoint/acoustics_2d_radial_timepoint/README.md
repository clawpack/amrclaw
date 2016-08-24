### Folder Organization
* **adjoint:**

Contains code to solve the adjoint problem.

The output times specified in this directory should agree with those for the
forward code. 

### Running the Code

* Go to the folder **adjoint** and run in a terminal:

```
make new
make .plots
```

The code will produce two new folders: _output and _plots. 
The first one contains all the output files, while the latter one contains the plots and interactive 
visualization apps.

* Go to the main folder **acoustics_2d_radial_timepoint** and run in the terminal:

```
make new
make .plots
```

The code will produce two new folders: _output and _plots. 
The first one contains all the output files, while the latter one contains the plots and interactive 
visualization apps.

### Running Variations

* Running the example with adjoint flagging:

Run in the terminal:

```
python run_adjoint_flagging.py
```

This example can be run with either the regular pressure-flagging technique (by using the flag2refine file) or with error-flagging (by using errf1 file). Either option can be set by setting flag_richardson and/or flag_richardson in setrun.py.
