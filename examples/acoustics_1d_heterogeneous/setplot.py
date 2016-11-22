def setplot(plotdata):
    
    plotdata.clearfigures()

    # Tuples of (variable name, variable number)
    figdata = [('Pressure', 0),
               ('Velocity', 1)]

    # Draw a vertical dashed line at the interface
    # between different media
    def draw_interface(current_data):
        import pylab
        pylab.plot([0., 0.], [-1000., 1000.], 'k--')
    
    # Formatting title
    def format(current_data, var, adjoint):
        import pylab
        size = 28
        t = current_data.t
        timestr = float(t)
        
        if (var == 0):
            titlestr = 'Pressure at time %s' % timestr
        else:
            titlestr = 'Velocity at time %s' % timestr
        
        pylab.title(titlestr, fontsize=size)
        pylab.xticks(fontsize=size)
        pylab.yticks(fontsize=size)
        pylab.xticks([-4, -2, 0, 2], fontsize=size)
    
    # After axis function for pressure
    def aa_pressure(current_data):
        draw_interface(current_data)
        format(current_data, 0, False)

    # After axis function for velocity
    def aa_velocity(current_data):
        draw_interface(current_data)
        format(current_data, 1, False)

    for varname, varid in figdata:
        plotfigure = plotdata.new_plotfigure(name=varname, figno=varid)
        plotfigure.kwargs = {'figsize': (11,5)}

        plotaxes = plotfigure.new_plotaxes()
        plotaxes.axescmd = 'axes([0.1,0.1,0.4,0.8])'
        plotaxes.xlimits = [-5., 3.]
        plotaxes.ylimits = [-0.5, 1.5]    # Good for both vars because of near-unit impedance
        if (varid == 0):
            plotaxes.afteraxes = aa_pressure
        else:
            plotaxes.afteraxes = aa_velocity

        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.plot_var = varid
        plotitem.color = 'b'
        plotitem.plotstyle = 'o'

    plotdata.printfigs = True          # Whether to output figures
    plotdata.print_format = 'png'      # What type of output format
    plotdata.print_framenos = 'all'    # Which frames to output
    plotdata.print_fignos = 'all'      # Which figures to print
    plotdata.html = True               # Whether to create HTML files
    plotdata.latex = False             # Whether to make LaTeX output

    return plotdata
