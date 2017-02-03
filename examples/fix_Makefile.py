"""
Fix Makefile for changes to support dimensional splitting, based
on amrclaw@260e3598e.

Modified for 3d
"""

from __future__ import absolute_import
from __future__ import print_function
import os

f = open('Makefile').read()
f = f.replace('stepgrid.f \\\n ', 'stepgrid.f \\\n  $(AMRLIB)/stepgrid_dimSplit.f \\\n ')
print(f.find('$(AMRLIB)/step3.f \\\n '))
f = f.replace('$(AMRLIB)/step3.f \\\n ', '$(AMRLIB)/step3.f \\\n  $(AMRLIB)/step3x.f \\\n  $(AMRLIB)/step3y.f \\\n  $(AMRLIB)/step3z.f \\\n ') 
f = f.replace('flux3.f \\\n ', 'flux3.f \\\n  $(AMRLIB)/flux3_dimSplit.f \\\n ')

os.system('mv Makefile Makefile_orig')
open('Makefile','w').write(f)
print('Modified Makefile, original in Makefile_orig')

