
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import absolute_import
from __future__ import print_function

from clawpack.clawutil.data import ClawData
from numpy import linspace
probdata = ClawData()
probdata.read('setprob.data', force=True)
print("Parameters: u = %g, beta = %g" % (probdata.u, probdata.beta))

def qtrue(x,t):
    """
    The true solution, for comparison.  
    Should be consistent with the initial data specified in qinit.f90.
    """
    from numpy import mod, exp, where, logical_and
    x0 = x - probdata.u*t
    x0 = mod(x0, 1.)   # because of periodic boundary conditions
    q = exp(-probdata.beta * (x0-0.75)**2)
    q = where(logical_and(x0 > 0.1, x0 < 0.4), q+1, q)
    return q
    

#--------------------------
def setplot(plotdata=None):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Solution', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [-.6,1.2]
    plotaxes.title = 'q'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.amr_color = ['g','b','r']
    plotitem.amr_plotstyle = ['^-','s-','o-']
    plotitem.amr_data_show = [1,1,1]
    plotitem.amr_kwargs = [{'markersize':8},{'markersize':6},{'markersize':5}]

    # Plot true solution for comparison:
    def plot_qtrue(current_data):
        from pylab import plot, legend
        x = linspace(0,1,1000)
        t = current_data.t
        q = qtrue(x,t)
        plot(x,q,'k',label='true solution')
    
    def plot_qtrue_with_legend(current_data):
        from pylab import plot, legend
        x = linspace(0,1,1000)
        t = current_data.t
        q = qtrue(x,t)
        plot(x,q,'k',label='true solution')
        try:
            from clawpack.visclaw import legend_tools
            labels = ['Level 1','Level 2', 'Level 3','True solution']
            legend_tools.add_legend(labels, colors=['g','b','r','k'],
                        markers=['^','s','o',''], linestyles=['','','','-'],
                        loc='lower right')
        except:
            legend(loc='lower right')

    plotaxes.afteraxes = plot_qtrue_with_legend

    # ------------------------------------------
    # Figure with each level plotted separately:

    plotfigure = plotdata.new_plotfigure(name='By AMR Level', figno=2)
    plotfigure.kwargs = {'figsize':(8,10)}


    for level in range(1,4):
        # Set up plot for this level:
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.axescmd = 'subplot(3,1,%i)' % level
        plotaxes.xlimits = [0,1]
        plotaxes.ylimits = [-.5,1.3]
        plotaxes.title = 'Level %s' % level
        plotaxes.afteraxes = plot_qtrue

        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.plot_var = 0
        plotitem.amr_color = ['g','b','r']
        plotitem.amr_plotstyle = ['^-','s-','o-']
        plotitem.amr_data_show = [0,0,0]
        plotitem.amr_data_show[level-1] = 1  # show only one level



    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                     type='each_gauge')
    plotfigure.clf_each_gauge = True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Solution'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
