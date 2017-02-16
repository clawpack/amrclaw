import numpy
import matplotlib
import matplotlib.pyplot as plt
import slice_tools
from clawpack.visclaw.data import ClawPlotData
from clawpack.visclaw.frametools import plotframe
import setplot

plotdata = ClawPlotData()
plotdata = setplot.setplot(plotdata)
plotdata.outdir = '_output'

for i in range(1,21):

    frameno = i
    plotframe(frameno, plotdata)

    frameno = i
    framesoln = plotdata.getframe(frameno)

    xpts = numpy.linspace(-1,1,51)
    qpts, level_pts = slice_tools.interp_amr_transect(framesoln, xpts)

    plt.figure(figsize=(8,8))
    plt.subplot(211)
    plt.plot(xpts,qpts,'bo')
    plt.ylim(-1,1.)
    plt.title("interpolated solution")

    plt.subplot(212)
    plt.plot(xpts,level_pts,'bo')
    plt.ylim(-0.5, 3)
    plt.title("AMR level")
    plt.show()
