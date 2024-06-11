"""
    Running adjoint example
"""

import os

current_dir = os.getcwd()
adjoint_dir = os.path.join(current_dir, "adjoint")

# Running the adjoint problem
os.chdir(adjoint_dir)
os.system('make new')
os.system('make .plots')

# Running the forward problem
os.chdir(current_dir)
os.system('make new')
os.system('make .plots')

print('Finished running example with adjoint refinement')