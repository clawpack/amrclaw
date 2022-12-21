"""
Tools for defining and working with ruled rectangles.
"""

import numpy as np

class RuledRectangle(object):
    
    def __init__(self, fname=None, slu=None, rect=None):
        self.ixy = None  # 1 or 'x' if s=x,   2 or 'y' if s=y
        self.s = None    # array
        self.lower = None  # array of same length as s
        self.upper = None  # array of same length as s
        self.method = 0  # 0 for pw constant, 1 for pw linear
        self.ds = -1   # > 0 if s values are equally spaced
        if fname is not None:
            self.read(fname)
            if slu is not None:
                print('*** Warning: ignoring slu since fname also given')
            if rect is not None:
                print('*** Warning: ignoring rect since fname also given')
        elif slu is not None:
            self.s = slu[:,0]
            self.lower = slu[:,1]
            self.upper = slu[:,2]
            if rect is not None:
                print('*** Warning: ignoring rect since slu also given')
            assert np.diff(self.s).min() > 0, \
                '*** s must be monotonically increasing: \n  s =  %s' % self.s
            assert np.all(slu[:,1] <= slu[:,2]), \
                '*** values in L column of slu array must be <= U column'
        elif rect is not None:
            # define a simple rectangle:
            x1,x2,y1,y2 = rect
            self.s = np.array([x1,x2])
            self.lower = np.array([y1,y1])
            self.upper = np.array([y2,y2])
            self.ixy = 1
            self.method = 0
            self.ds = x2 - x1
        
    def bounding_box(self):
        if self.ixy in [1,'x']:
            x1 = self.s.min()
            x2 = self.s.max()
            y1 = self.lower.min()
            y2 = self.upper.max()
        elif self.ixy in [2,'y']:
            y1 = self.s.min()
            y2 = self.s.max()
            x1 = self.lower.min()
            x2 = self.upper.max()
        else:
            raise ValueError('Unrecognized attribute ixy = %s' % self.ixy)
        return [x1, x2, y1, y2]

    def slu(self):
        """
        Return npts by 3 array with columns (s, lower, upper)
        """
        assert np.diff(self.s).min() > 0, \
            '*** s must be monotonically increasing: \n  s =  %s' % self.s
        slu = np.vstack((self.s, self.lower, self.upper)).T
        return slu

            
    def vertices(self):
        """
        Return the vertices (x,y) of the polygon defined by this RuledRectangle.
        """
        
        assert np.diff(self.s).min() > 0, \
            '*** s must be monotonically increasing: \n  s =  %s' % self.s
        assert np.all(self.lower <= self.upper), \
            '*** values in lower array should be <= values in upper array'
            
        if self.method == 0:
            # stacked rectangles, requires doubling up points for vertices
            ss = [self.s[0]]
            for sk in self.s[1:-1]:
                ss += [sk,sk]
            ss = np.array(ss + [self.s[-1]])
            ss = np.hstack((ss, np.flipud(ss), self.s[0]))
            ll = []
            uu = []
            for k in range(len(self.s)-1):
                ll += [self.lower[k],self.lower[k]]
                uu += [self.upper[k],self.upper[k]]
            lu = np.hstack((ll, np.flipud(uu), self.lower[0]))
            
        elif self.method == 1:
            # points defined by s, lower, upper are connected by lines
            ss = np.hstack((self.s, np.flipud(self.s), self.s[0]))
            lu = np.hstack((self.lower, np.flipud(self.upper), self.lower[0]))
            
        if self.ixy in [1,'x']:
            x = ss
            y = lu
        elif self.ixy in [2,'y']:
            x = lu
            y = ss
        else:
            raise ValueError('Unrecognized attribute ixy = %s' % self.ixy)
        return x,y
            
            
    def mask_outside(self, X, Y):
        """
        Given 2d arrays X,Y, return a mask with the same shape with
        mask == True at points that are outside this RuledRectangle.
        So if Z is a data array defined at points X,Y then 
            ma.masked_array(Z, mask) 
        will be a masked array that can be used to plot only the values
        inside the Ruled Region.
        
        Works for self.method == 0 or 1, and also
              for self.s monotonically increasing or decreasing, 
                         but elsewhere we require increasing so check for that.
        
        """
        assert np.diff(self.s).min() > 0, \
            '*** s must be monotonically increasing: \n  s =  %s' % self.s
        assert np.all(self.lower <= self.upper), \
            '*** values in lower array should be <= values in upper array'
        
        transpose_arrays =  (X[0,0] == X[0,-1])
        if transpose_arrays:
            x = X[:,0]
            y = Y[0,:]
        else:
            x = X[0,:]
            y = Y[:,0]
        assert x[0] != x[-1], '*** Wrong orientation?'            
        
        def bracket(xy, s):
            """
            Given a point xy (an x or y value)
            and an array s that is monotone increasing or decreasing, 
            return indices (k1, k2) such that 
                s[k1] <= xy <= s[k2]
            If xy is outside the range of s, return (-1,-1)
            """
            from numpy import ma, argmin
            if xy < s.min() or xy > s.max():
                return -1,-1
            pos = ma.masked_where(xy-s < 0, xy-s)
            k1 = argmin(pos)
            pos = ma.masked_where(s-xy < 0, s-xy)
            k2 = argmin(pos)
            return k1,k2
            
        mask = np.empty((len(y),len(x)), dtype=bool)
        mask[:,:] = True
        if self.ixy in [1,'x']:
            # indices of x array that lie within s.min() to s.max():
            iin, = np.where(np.logical_and(self.s.min() <= x, x <= self.s.max()))
            
            for i in iin:
                # a column of array that might include points in the polygon
                xi = x[i]
                # determine k1,k2 such that s[k1] <= xi <= s[k2]:
                k1,k2 = bracket(xi, self.s)
                
                if k1 > -1:
                    sk1 = self.s[k1]
                    sk2 = self.s[k2]
                    if (self.method==0) or (sk1 >= xi):
                        ylower = self.lower[k1]
                        yupper = self.upper[k1]
                    else:
                        # method==1 and not at top, so linear interpolation:
                        alpha = (xi-sk1)/(sk2-sk1)
                        ylower = (1-alpha)*self.lower[k1] + \
                                     alpha*self.lower[k2]
                        yupper = (1-alpha)*self.upper[k1] + \
                                     alpha*self.upper[k2] 
                    # indices j for which ylower <= y[j] <= yupper, so inside:  
                    j, = np.where(np.logical_and(ylower <= y, y <= yupper))
                    mask[j,i] = False
                    
        elif self.ixy in [2,'y']:
            # indices of y array that lie within s.min() to s.max():
            jin, = np.where(np.logical_and(self.s.min() <= y, y <= self.s.max()))
            
            for j in jin:
                # a row of array that might include points in the polygon
                yj = y[j]
                # determine k1,k2 such that s[k1] <= yj <= s[k2]:
                k1,k2 = bracket(yj, self.s)
                
                if k1 > -1:
                    sk1 = self.s[k1]
                    sk2 = self.s[k2]
                    if (self.method==0) or (min(sk1,sk2) >= yj):
                        xlower = self.lower[k1]
                        xupper = self.upper[k1]
                    else:
                        # method==1 and not at top, so linear interpolation:
                        alpha = (yj-sk1)/(sk2-sk1)
                        xlower = (1-alpha)*self.lower[k1] + \
                                     alpha*self.lower[k2]
                        xupper = (1-alpha)*self.upper[k1] + \
                                     alpha*self.upper[k2]   
                    # indices i for which xlower <= x[i] <= xupper, so inside:
                    i, = np.where(np.logical_and(xlower <= x, x <= xupper))
                    mask[j,i] = False

        else:
            raise ValueError('Unrecognized attribute ixy = %s' % self.ixy)
                    
        if transpose_arrays:
            mask = mask.T
            
        return mask

                    
    def write(self, fname, verbose=False):
        slu = self.slu()
        ds = self.s[1:] - self.s[:-1]
        dss = ds.max() - ds.min()
        if dss < 1e-6*ds.max():
            self.ds = ds.max()
        else:
            self.ds = -1  # not uniformly spaced
            
        # if ixy is 'x' or 'y' replace by 1 or 2 for writing:
        if self.ixy in [1,'x']:
            ixyint = 1
        elif self.ixy in [2,'y']:
            ixyint = 2
        else:
            raise ValueError('Unrecognized attribute ixy = %s' % self.ixy)
            
        header = """\n%i   ixy\n%i   method\n%g    ds\n%i    nrules""" \
            % (ixyint,self.method,self.ds,len(self.s))
        np.savetxt(fname, slu,header=header,comments='',fmt='%.9f  ')
        
        if verbose:
            print("Created %s" % fname)

    def read(self, fname):
        lines = open(fname,'r').readlines()
        k = -1
        comments = True
        while comments:
            k += 1
            line = lines[k].strip()
            if (line != '') and (line[0] != '#'):
                comments = False
        self.ixy = int(line.split()[0])
        k += 1
        self.method = int(lines[k].split()[0])
        k += 1
        self.ds = float(lines[k].split()[0])
        k += 1
        self.nrules = int(lines[k].split()[0])
        slu = np.loadtxt(fname, skiprows=k+1)
        assert slu.shape[0] == self.nrules, '*** wrong shape'
        self.s = slu[:,0]
        self.lower = slu[:,1]
        self.upper = slu[:,2]
        assert np.all(slu[:,1] <= slu[:,2]), \
            '*** values in L column of slu array must be <= U column'
        
    def make_kml(self, fname='RuledRectangle.kml', name='RuledRectangle', 
                 color='00FFFF', width=2, verbose=False):
        from clawpack.geoclaw import kmltools
        x,y = self.vertices()
        kmltools.poly2kml((x,y), fname=fname, name=name, color=color, 
                          width=width, verbose=verbose)
        
        
def ruledrectangle_covering_selected_points(X, Y, pts_chosen, ixy, method=0,
                                            padding=0, verbose=True):
    """
    Given 2d arrays X,Y and an array pts_chosen of the same shape,
    returns a RuledRectangle that covers these points as compactly as possible.
    If method==0 then the RuledRectangle will be "stacked rectangles" that 
    cover all the cells indicated by pts_chosen 
    (assuming the X,Y points are cell centers on a uniform grid).
    
    If method==1, then the RuledRectangle will cover the cell centers but
    not all the full grid cells, since the edges will be lines connecting
    some cell centers.
    """
        
    if padding != 0:
        print('*** Warning, padding != 0 is not implemented!')
        
    if np.ndim(X) == 2:
        x = X[0,:]
        y = Y[:,0]
    else:
        x = X
        y = Y

    dx = x[1] - x[0]
    dy = y[1] - y[0]

    if ixy in [1,'x']:

        # Ruled rectangle with s = x:

        s = []
        lower = []
        upper = []
        for i in range(len(x)):
            if pts_chosen[:,i].sum() > 0:
                j = np.where(pts_chosen[:,i]==1)[0]
                j1 = j.min()
                j2 = j.max()
                s.append(x[i])
                lower.append(y[j1])
                upper.append(y[j2])
                
    elif ixy in [2,'y']:

        # Ruled rectangle with s = y:

        s = []
        lower = []
        upper = []

        for j in range(len(y)):
            if pts_chosen[j,:].sum() > 0:
                i = np.where(pts_chosen[j,:]==1)[0]
                i1 = i.min()
                i2 = i.max()
                s.append(y[j])
                lower.append(x[i1])
                upper.append(x[i2])
                
    else:
        raise(ValueError('Unrecognized value of ixy'))

    s = np.array(s)
    lower = np.array(lower)
    upper = np.array(upper)
    ds = s[1] - s[0]

    if method == 0:
        # extend so rectangles cover grid cells with centers at (x,y)
        if verbose:
            print('Extending rectangles to cover grid cells')
        if abs(dx - dy) > 1e-6:
            print('*** Warning, dx = %.8e not equal to dy = %.8e' \
                  % (dx, dy))
        lower = lower - 0.5*ds
        upper = upper + 0.5*ds
        s = s - 0.5*ds
        s = np.hstack((s, s[-1]+ds))
        lower = np.hstack((lower, lower[-1]))
        upper = np.hstack((upper, upper[-1]))
        
    rr = RuledRectangle()
    rr.ixy = ixy
    rr.s = s
    rr.lower = lower
    rr.upper = upper
    rr.method = method
    rr.ds = ds

    if verbose:
        rr_npts = int(np.ceil(np.sum(rr.upper - rr.lower) / ds))
        print('RuledRectangle covers %s grid points' % rr_npts)

    return rr
