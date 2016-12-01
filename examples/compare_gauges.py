"""
Compare gauges between two different runs, useful for regression testing.

Requires visclaw@c3bc5aab82 or later for gaugetools.compare_gauges.
"""

from __future__ import absolute_import
from __future__ import print_function
from clawpack.visclaw import gaugetools

outdir1 = '_output_original'
outdir2 = '_output'

tol = 0.  # tolerance for comparison
    
matches = gaugetools.compare_gauges(outdir1, outdir2, 
            gaugenos='all', q_components='all',
            tol=tol, verbose=True, plot=False)

if not matches:
    print("*** Warning: results to not all match to tol = %g" % tol)
    print("*** You might want to call gaugetools.compare_gauges with ")
    print("    plot=True to view differences.")


