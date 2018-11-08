Acoustics 1D -- Heterogeneous example
------------------------------------------

1D acoustics in a piecewise constant medium to illustrate reflection and
transmission at an interface.

The density and bulk modulus of the medium are specified in setrun.py, along with a parameter determining what initial conditions to use (see qinit.f).

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
The first one contains all the output files, while the latter one contains the plots and interactive visualization apps.

* Go to the main folder **acoustics_1d_adjoint** and run in the terminal:

```
make new
make .plots
```

The code will produce two new folders: _output and _plots. 
The first one contains all the output files, while the latter one contains the plots and interactive visualization apps.

### Running Variations

* Running the example with adjoint flagging:

Run in the terminal:

```
python run_adjoint_flagging.py
```

This example can be run with either the regular surface-flagging technique (by using the flag2refine file) or with error-flagging (by using errf1 file). These options can be set by setting the flag_forward_adjoint and/or flag_richardson_adjoint flags in setrun.py.
