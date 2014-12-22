"""
Fix Makefile for changes to support dimensional splitting, based
on amrclaw@260e3598e.
"""

import os

f = open('Makefile').read()
f = f.replace('stepgrid.f \\\n ', 'stepgrid.f \\\n  $(AMRLIB)/stepgrid_dimSplit.f \\\n ')
print f.find('$(AMRLIB)/step2.f90 \\\n ')
f = f.replace('$(AMRLIB)/step2.f90 \\\n ', '$(AMRLIB)/step2.f90 \\\n  $(AMRLIB)/step2x.f90 \\\n  $(AMRLIB)/step2y.f90 \\\n ')
f = f.replace('flux2.f \\\n ', 'flux2.f \\\n  $(AMRLIB)/flux2_dimSplit.f \\\n ')

os.system('mv Makefile Makefile_orig')
open('Makefile','w').write(f)
print 'Modified Makefile, original in Makefile_orig'

