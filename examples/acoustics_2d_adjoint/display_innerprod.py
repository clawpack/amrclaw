from pylab import *
from clawpack.visclaw.data import ClawPlotData
from clawpack.visclaw.frametools import plotframe
from clawpack.visclaw import colormaps
import setplot

def main():

    # Change these values, and you should not have to
    # modify the rest of the file
    # ----------------------------------
    
    # Frames that you wish to consider
    # (the inner product will be taken between these two frames)
    forward_frameno = 8
    adjoint_frameno = 7
    
    # Number of levels in the simulation
    levels = [1,2,3]
    # ----------------------------------

    plotdata_f = ClawPlotData()
    plotdata_f = setplot.setplot(plotdata_f)
    plotdata_f.outdir = '_output'

    plotdata_a = ClawPlotData()
    plotdata_a = setplot.setplot(plotdata_a)
    plotdata_a.outdir = 'adjoint/_output'

    # Read in frame from forward solution
    framesoln_f = plotdata_f.getframe(forward_frameno)

    # Read in frame from adjoint solution
    framesoln_a = plotdata_a.getframe(adjoint_frameno)

    white_red = colormaps.make_colormap({0.0:'w', 0.2:'r'})

    states_f = framesoln_f.states

    for current_level in levels:
        # loop over all patches in forward solution
        for state_f in states_f:
            patch = state_f.patch
            level = patch.level
            if (level == current_level):
                q = state_f.q
                xlower = patch.dimensions[0].lower
                xupper = patch.dimensions[0].upper
                ylower = patch.dimensions[1].lower
                yupper = patch.dimensions[1].upper
                dx = patch.delta[0]
                dy = patch.delta[1]
                meqn,mx,my = q.shape
                xc_edges, yc_edges = patch.grid.c_nodes
                xc_centers, yc_centers = patch.grid.c_centers
    
                # Creating array to hold inner product
                innerprod = np.zeros(q[2,:,:].shape)
    
                # loop over each cell in this patch of the forward solution
                for (i,row) in enumerate(q[0,:,:]):
                    for (j,val) in enumerate(row):
                        x = xc_centers[i,j]
                        y = yc_centers[i,j]

                        calculate_innerprod(q,innerprod,x,y,i,j,framesoln_a)
                    # done looping over cells in row
                # done looping over rows in patch
                pcolor(xc_edges, yc_edges, innerprod, cmap=white_red)
            # done with patches on current level
        # done looping over patches
    # done looping over levels

    axis([-4, 8, -1, 11])
    show()

#------------------------------------------------------------------
def calculate_innerprod(q,innerprod,x,y,i,j,framesoln_a):
#------------------------------------------------------------------

    # loop over all patches in adjoint solution
    states_a = framesoln_a.states
    for state_a in states_a:
        q_a = state_a.q
        patch_a = state_a.patch
        xlower_a = patch_a.dimensions[0].lower
        xupper_a = patch_a.dimensions[0].upper
        ylower_a = patch_a.dimensions[1].lower
        yupper_a = patch_a.dimensions[1].upper
        dx_a = patch_a.delta[0]
        meqn_a,mx_a,my_a = q_a.shape
        xc_edges_a, yc_edges_a = patch_a.grid.c_nodes
        xc_centers_a, yc_centers_a = patch_a.grid.c_centers


    if ((x > xlower_a) and (x < xupper_a)
        and (y > ylower_a) and (y < yupper_a)):
        # finding correct cells to interpolate with
        idx = (np.abs(xc_centers_a[:,0]-x)).argmin()
        idy = (np.abs(yc_centers_a[0,:]-y)).argmin()
            
        if (xc_centers_a[idx,0] <= x):
            idxm = idx
            idxp = idx + 1
            if (idxp >= mx_a):
                idxp = idx
        else:
            idxm = idx - 1
            idxp = idx
            if (idxm < 0):
                idxm = idx
                
        if (yc_centers_a[0,idy] <= y):
            idym = idy
            idyp = idy + 1
            if (idyp >= my_a):
                idyp = idy
        else:
            idym = idy - 1
            idyp = idy
            if (idym < 0):
                idym = idy
                
        # interpolate in x
        xinterp = False
        if (idxm != idxp):
            xinterp = True
                        
            xp = xc_centers_a[idxp,0]
            xm = xc_centers_a[idxm,0]
            denomx = xp - xm
                        
            q_interp1 = ((xp - x)/denomx)*q_a[:,idxp,idyp] + ((x - xm)/denomx)*q_a[:,idxm,idyp]
            q_interp2 = ((xp - x)/denomx)*q_a[:,idxp,idym] + ((x - xm)/denomx)*q_a[:,idxm,idym]
                
        # interpolate in y
        if (idym != idyp):
            yp = yc_centers_a[0,idyp]
            ym = yc_centers_a[0,idym]
            denomy = yp - ym
                        
            if (xinterp):
                q_interp = ((yp - y)/denomy)*q_interp1 + ((y - ym)/denomy)*q_interp2
            else:
                q_interp = ((yp - y)/denomy)*q_a[:,idxp,idyp] + ((y - ym)/denomy)*q_a[:,idxp,idym]
        else:
            if (xinterp):
                q_interp = q_interp1
            else:
                q_interp = q_a[:,idxp,idyp]
                            
        # calculating inner product
                            
        # leaving off the last term in q, because that is the inner product
        # that was calculated during the actual Clawpack run
        innerprod[i,j] = abs(np.sum(q[:-1,i,j]*q_interp))

if __name__ == '__main__':
    main()
