
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

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

    gamma = 1.4
    def pressure(current_data):
        q = current_data.q
        rho = q[0,:]
        u = q[1,:]/rho
        p = (gamma-1)*(q[2,:] - 0.5*rho*u**2)
        return p

    level_labels = ['Level 1','Level 2', 'Level 3', 'Level 4']
    level_colors = ['k','g','b','r']
    level_styles = '-'
    level_markers = ''

    def add_legend(current_data):
        try:
            from clawpack.visclaw import legend_tools
            legend_tools.add_legend(labels=level_labels, colors=level_colors,
                        markers=level_markers, linestyles=level_styles,
                        loc='upper right')
        except:
            print('legend_tools error')
            pass


    # Density and pressure plots
    # --------------------------

    plotfigure = plotdata.new_plotfigure(name='Density and Pressure', figno=1)
    plotfigure.kwargs = {'figsize':(10,5)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(1,2,1)'   # left figure
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,7]
    plotaxes.title = 'Density'
    #plotaxes.afteraxes = add_legend

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_pwconst')
    plotitem.plot_var = 0
    plotitem.amr_color = level_colors
    plotitem.amr_plotstyle = level_styles
    plotitem.amr_data_show = [1]  # all levels


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(1,2,2)'   # right
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [-10,1100]
    plotaxes.title = 'Pressure'
    plotaxes.afteraxes = add_legend

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_pwconst')
    plotitem.plot_var = pressure
    plotitem.amr_color = level_colors
    plotitem.amr_plotstyle = level_styles
    plotitem.amr_data_show = [1]  # all levels

    #  Density and pressure plots, zoomed at final time
    # -------------------------------------------------

    plotfigure = plotdata.new_plotfigure(name='zoom view', figno=2)
    plotfigure.kwargs = {'figsize':(10,5)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(1,2,1)'   # left figure
    plotaxes.xlimits = [0.5,0.9]
    plotaxes.ylimits = [0,7]
    plotaxes.title = 'Density'
    #plotaxes.afteraxes = add_legend

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_pwconst')
    plotitem.plot_var = 0
    plotitem.amr_color = level_colors
    plotitem.amr_plotstyle = level_styles
    plotitem.amr_data_show = [1]  # all levels


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(1,2,2)'   # right
    plotaxes.xlimits = [0.5,0.9]
    plotaxes.ylimits = [-10,450]
    plotaxes.title = 'Pressure'
    plotaxes.afteraxes = add_legend

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_pwconst')
    plotitem.plot_var = pressure
    plotitem.amr_color = level_colors
    plotitem.amr_plotstyle = level_styles
    plotitem.amr_data_show = [1]  # all levels

    
    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                         type='each_gauge')
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-1,30]
    plotaxes.title = 'Density'
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

    
