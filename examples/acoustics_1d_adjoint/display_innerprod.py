
from pylab import *
from clawpack.visclaw.data import ClawPlotData
from clawpack.visclaw.frametools import plotframe
import setplot
import matplotlib.lines as mlines

forward_frameno = 5
adjoint_frameno = 16

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

states_f = framesoln_f.states
# loop over all patches in forward solution
for state_f in states_f:
    q = state_f.q
    patch = state_f.patch
    xlower = patch.dimensions[0].lower
    xupper = patch.dimensions[0].upper
    dx = patch.delta[0]
    meqn,mx = q.shape
    x_vals = state_f.grid.x.centers
    level = patch.level
    
    # Creating array to hold inner product
    innerprod = np.zeros(q[2,:].shape)

    # loop over all patches in adjoint solution
    states_a = framesoln_a.states
    for state_a in states_a:
        q_a = state_a.q
        patch_a = state_a.patch
        xlower_a = patch_a.dimensions[0].lower
        xupper_a = patch_a.dimensions[0].upper
        dx_a = patch_a.delta[0]
        meqn_a,mx_a = q_a.shape
        x_vals_a = state_a.grid.x.centers
    
        # loop over each cell in this patch of the forward solution
        counter = 0
        for val in q[0,:]:
            x = x_vals[counter]
            if ((x > xlower_a) and (x < xupper_a)):
                # finding correct cells to interpolate with
                idx = (np.abs(x_vals_a-x)).argmin()

                if (x_vals_a[idx] <= x):
                    idxm = idx
                    idxp = idx + 1
                elif (x_vals_a[idx] > x):
                    idxm = idx - 1
                    idxp = idx
        
                # interpolate between these two cells
                xp = x_vals_a[idxp]
                xm = x_vals_a[idxm]
                denom = xp - xm
                    
                q_interp = ((xp - x)/denom)*q_a[:,idxp] + ((x - xm)/denom)*q_a[:,idxm]
                
                # calculating inner product
                innerprod[counter] = q[0,counter] * q_interp[0] + q[1,counter]*q_interp[1]
            counter = counter + 1

    # Plotting the inner product for this frame
    if level == 1:
        plot(x_vals, innerprod, '^-g')
    if level == 2:
        plot(x_vals, innerprod, 's-b')

labels = ['Level 1','Level 2']
legend(labels, loc='upper left')
axis([-5, 3, -.5, 1.1])
plot([0., 0.], [-1000., 1000.], 'k--')
show()
