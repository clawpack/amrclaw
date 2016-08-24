8
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
"""
#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """
    from clawpack.visclaw import colormaps
    
    # Reversing time in adjoint output
    setadjoint()

    plotdata.clearfigures()  # clear any old figures,axes,items data
    

    # Figure for pressure
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=0)
    plotfigure.kwargs = {'figsize': (14,5)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([-0.22,0.1,0.8,0.8])'
    plotaxes.xlimits = [-3,8]
    plotaxes.ylimits = [-1,10]
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = fixup

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.add_colorbar = False
    plotitem.show = True       # show on plot?
    plotitem.pcolor_cmin = -1.0
    plotitem.pcolor_cmax = 1.0
    plotitem.amr_patchedges_show = [1,1,1]
    plotitem.amr_celledges_show = [0,0,0]
    
    # Adding innerproduct plot
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Inner Product'
    plotaxes.axescmd = 'axes([0.11,0.1,0.8,0.8])'
    plotaxes.xlimits = [-3,8]
    plotaxes.ylimits = [-1,10]
    plotaxes.title = 'Inner Product'
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = fixup_innerprod
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 3
    plotitem.pcolor_cmap = colormaps.white_red
    plotitem.add_colorbar = False
    plotitem.show = True       # show on plot?
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 0.15
    plotitem.amr_patchedges_show = [0,0,0]
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_data_show = [1,1,0]
    
    # Adding adjoint plot
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Inner Product'
    plotaxes.axescmd = 'axes([0.44,0.1,0.8,0.8])'
    plotaxes.xlimits = [-3,8]
    plotaxes.ylimits = [-1,10]
    plotaxes.title = 'Adjoint'
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = fixup_adjoint
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    import os
    plotitem.outdir = os.path.join(os.getcwd(), 'adjoint/_outputReversed')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.add_colorbar = False
    plotitem.show = True       # show on plot?
    plotitem.pcolor_cmin = -0.2
    plotitem.pcolor_cmax = 0.2
    plotitem.amr_patchedges_show = [0,0,0]
    plotitem.amr_celledges_show = [0,0,0]
    
    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True
                                         
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Pressure'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'
    plotitem.kwargs = {'linewidth': 3}
    plotaxes.afteraxes = fixup_gauge

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

# Afteraxis functions:

def addgauges(current_data):
    from clawpack.visclaw import gaugetools
    gaugetools.plot_gauge_locations(current_data.plotdata, \
         gaugenos='all', format_string='ko', add_labels=True)

def fixup(current_data):
    import pylab
    size = 34
    addgauges(current_data)
    pylab.title('Forward Pressure', fontsize=size)
    pylab.xticks([-2, 0, 2, 4, 6], fontsize=size)
    pylab.yticks([0, 2, 4, 6, 8], fontsize=size)

def fixup_innerprod(current_data):
    import pylab
    size = 36
    addgauges(current_data)
    pylab.title('Inner Product', fontsize=size)
    pylab.xticks([-2, 0, 2, 4, 6], fontsize=size)
    pylab.tick_params(axis='y', labelleft='off')

def fixup_adjoint(current_data):
    import pylab
    size = 36
    addgauges(current_data)
    pylab.title('Adjoint Pressure', fontsize=size)
    pylab.xticks([-2, 0, 2, 4, 6], fontsize=size)
    pylab.tick_params(axis='y', labelleft='off')

def fixup_gauge(current_data):
    import pylab
    size = 36
    pylab.title('Pressure at Gauge 0', fontsize=size)
    pylab.xticks([1.25, 1.35, 1.45], fontsize=size)
    pylab.yticks([0, 0.1, 0.2, 0.3, 0.4], fontsize=size)


#-------------------
def setadjoint():
    #-------------------
    
    """
        Reverse order of adjoint images, for plotting
        adjacent to forward plots.
        """
    
    import os,sys,glob
    from clawpack.pyclaw import io
    
    outdir = 'adjoint/_output'
    outdir2 = 'adjoint/_outputReversed'
    
    os.system('mkdir -p %s' % outdir2)
    
    files = glob.glob(outdir+'/fort.q0*')
    files.sort()
    n = len(files)
    
    if (n >= 1):
        # Find the final time.
        fname = files[n-1]
        fname = fname.replace('q','t')
        f = open(fname,'r')
        tfinal,meqn,npatches,maux,num_dim = io.ascii.read_t(n-1,path=outdir)
    
        for k in range(n):
            # Creating new files
            fname = files[k]
            newname = outdir2 + '/fort.q%s' % str(n-k-1).zfill(4)
            cmd = 'cp %s %s' % (fname,newname)
            os.system(cmd)
        
            fname = fname.replace('q','t')
            newname = newname.replace('q','t')
            cmd = 'cp %s %s' % (fname,newname)
            os.system(cmd)
        
            # Reversing time
            f = open(newname,'r+')
            frameno = n-k-1
            t,meqn,npatches,maux,num_dim = io.ascii.read_t(frameno,path=outdir2)
            t = tfinal - t
        
            # Writting new time out to file
            f.write('%18.8e     time\n' % t)
            f.write('%5i                  num_eqn\n' % meqn)
            f.write('%5i                  nstates\n' % npatches)
            f.write('%5i                  num_aux\n' % maux)
            f.write('%5i                  num_dim\n' % num_dim)
            f.close()
# end of function setadjoint
# ----------------------
