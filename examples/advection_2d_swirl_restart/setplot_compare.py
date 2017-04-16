
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import os
outdir1 = os.path.abspath('_output')
outdir2 = os.path.abspath('_output_restart')

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    

    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)
    plotfigure.kwargs={'figsize':(12,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(121)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'original'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.outdir = outdir1
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 1.0
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.5
    plotitem.amr_celledges_show = [True, True, False]
    plotitem.amr_patchedges_show = [True, True]
    plotitem.show = True       # show on plot?
    plotaxes = plotfigure.new_plotaxes()

    plotaxes.axescmd = 'subplot(122)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'restart'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.outdir = outdir2
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 1.0
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.5
    plotitem.amr_celledges_show = [True, True, False]
    plotitem.amr_patchedges_show = [True, True]
    plotitem.show = True       # show on plot?
    

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.,4]
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q'

    # Plot q as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.outdir = outdir1
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'

    # Plot restart q as red curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.outdir = outdir2
    plotitem.plot_var = 0
    plotitem.plotstyle = 'r-'

    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

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

    
