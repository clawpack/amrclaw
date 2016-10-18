
def interp_amr_transect(framesoln, xpts):

    from numpy import ma
    import numpy as np

    q_pts = ma.masked_array(np.zeros(xpts.shape),mask=True)
    level_pts = ma.masked_array(np.zeros(xpts.shape),mask=True)

    states = framesoln.states
    # loop over all patches
    for state in states:
        q = state.q
        patch = state.patch
        xlower = patch.dimensions[0].lower
        xupper = patch.dimensions[0].upper
        dx = patch.delta[0]
        meqn,mx = q.shape
        pts_in_patch = np.logical_and(xpts >= xlower, xpts <= xupper)

        level_pts[pts_in_patch] = patch.level

        q_in_patch = []
        for xp in xpts[pts_in_patch]:
            i = min(int(np.floor((xp - xlower)/dx)), mx-1)
            q_in_patch.append(q[0,i])
            
        q_pts[pts_in_patch] = q_in_patch

    return q_pts, level_pts
