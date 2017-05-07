
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
"""

from __future__ import absolute_import
from clawpack.clawutil.data import ClawData
probdata = ClawData()
probdata.read('setprob.data', force=True)
ubar = probdata.u
vbar = probdata.v
vbar = probdata.w

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps

    plotdata.clearfigures() # clear any old figures,axes,items data

    
    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,50]
    plotaxes.ylimits = [0,50]
    plotaxes.title = 'Solution'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = -150.
    plotitem.pcolor_cmax = 150.
    plotitem.add_colorbar = True

    plotitem.amr_celledges_show = [0,0]  
    plotitem.amr_patchedges_show = [0,1]

    
    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='error pcolor', figno=10)

    def error(current_data):
        x = current_data.x
        y = current_data.y
        z = current_data.z
        t = current_data.t
        q = current_data.q
        qtrue = (x-ubar*t) + (y-vbar*t) + (z-wbar*t)
        e = q[0,:,:,:] - qtrue
        #import pdb; pdb.set_trace()
        return e

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,50]
    plotaxes.ylimits = [0,50]
    plotaxes.zlimits = [0,50]
    plotaxes.title = 'Error'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = error
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = -0.000001
    plotitem.pcolor_cmax = 0.000001
    plotitem.add_colorbar = True

    plotitem.amr_celledges_show = [0,0]  
    plotitem.amr_patchedges_show = [0,1]


    # Figure for contour plot
    plotfigure = plotdata.new_plotfigure(name='contour', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,50]
    plotaxes.ylimits = [0,50]
    plotaxes.zlimits = [0,50]
    plotaxes.title = 'Solution'
    plotaxes.scaled = True
    #plotaxes.afteraxes = addgauges

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 0
    plotitem.contour_nlevels = 40
    plotitem.contour_min = -150.99
    plotitem.contour_max = 150.99
    plotitem.amr_contour_colors = ['r','g','b']  # color on each level
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 1 


    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,50]
    plotaxes.ylimits = [0,50]
    plotaxes.title = 'Grid patches'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1]

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.zlimits = 'auto'
    plotaxes.title = 'q'

    # Plot q as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'

    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.html_movie = 'JSAnimation'      # new style, or "4.x" for old style
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
    
# To plot gauge locations on pcolor or contour plot, use this as
# an afteraxis function:

def addgauges(current_data):
    from clawpack.visclaw import gaugetools
    gaugetools.plot_gauge_locations(current_data.plotdata, \
         gaugenos='all', format_string='ko', add_labels=True)
