"""
Fix Makefile for changes to gauges_module.f90: remove dumpgauges.f
on amrclaw@77612a61b3
"""

import os

f = open('Makefile').read()
f = f.replace('  $(AMRLIB)/dumpgauge.f \\\n','')

os.system('mv Makefile Makefile_orig')
open('Makefile','w').write(f)
print 'Modified Makefile, original in Makefile_orig'

