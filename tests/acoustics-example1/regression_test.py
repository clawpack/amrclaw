import os
from clawpack.clawutil.chardiff import chardiff_dir
from clawpack.clawutil.imagediff import imagediff_dir

outdir = '_output_new'
plotdir = '_plots_new'

regression_dir= "/Users/rjl/git/rjleveque/amrclaw_regression_13feb2013/acoustics-example1"

make_commands = """
    make new
    make .plots OUTDIR=%s PLOTDIR=%s
    """ % (outdir,plotdir)

if True:
    os.system(make_commands)

regression_outdir = regression_dir + '/_output'
regression_plotdir = regression_dir + '/_plots'

print "\n======= output comparison ======="
chardiff_dir(regression_outdir, outdir, dir3='_regression_output_cdiff')
print "\n======= plots comparison ======="
imagediff_dir(regression_plotdir, plotdir, dir3='_regression_plots_idiff')
