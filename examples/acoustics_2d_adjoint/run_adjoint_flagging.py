"""
    Running adjoint example
"""

import os

currentdir = os.getcwd()
adjointdir = currentdir + '/adjoint'

# Running the adjoint problem
os.chdir(adjointdir)
os.system('make new')
os.system('make .plots')

# Running the forward problem
os.chdir(currentdir)
os.system('make new')
os.system('make .plots')

print 'Finished running example with adjoint refinement'