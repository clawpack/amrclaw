
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import absolute_import
from __future__ import print_function

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

    def draw_interface_add_legend(current_data):
        from pylab import plot
        plot([0., 0.], [-1000., 1000.], 'k--')
        try:
            from clawpack.visclaw import legend_tools
            labels = ['Level 1','Level 2', 'Level 3']
            legend_tools.add_legend(labels, colors=['g','b','r'],
                        markers=['^','s','o'], linestyles=['','',''],
                        loc='upper left')
        except:
            pass


    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Pressure and Velocity', figno=1)
    plotfigure.kwargs = {'figsize': (8,8)}
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'   # top figure
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-.5,1.1]
    plotaxes.title = 'Pressure'
    plotaxes.afteraxes = draw_interface_add_legend

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.amr_color = ['g','b','r']
    plotitem.amr_plotstyle = ['^-','s-','o-']
    plotitem.amr_data_show = [1,1,1]
    plotitem.amr_kwargs = [{'markersize':5},{'markersize':4},{'markersize':3}]

    # Figure for q[1]

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'   # bottom figure
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-.5,1.1]
    plotaxes.title = 'Velocity'
    plotaxes.afteraxes = draw_interface_add_legend

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 1
    plotitem.amr_color = ['g','b','r']
    plotitem.amr_plotstyle = ['^-','s-','o-']
    plotitem.amr_data_show = [1,1,1]
    plotitem.amr_kwargs = [{'markersize':5},{'markersize':4},{'markersize':3}]

    
    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                         type='each_gauge')
    plotfigure.clf_each_gauge = True

                                         
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Pressure'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Velocity'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 1
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

    
