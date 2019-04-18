import numpy as np
import bisect

class Corner(object):
    
    def __init__(self, side, x, y, yy):
        self.side = side
        self.x = x
        self.y = y
        self.yy = yy

class Sweep(object):
    
    def __init__(self, xmin, ymin, ymax):
        self.x = xmin
        self.y = [ymin, ymax]
        self.n = [0]
        
    def update(self, x):
        self.x = x
        
    def split(self, corner):
        ind1 = (bisect.bisect_left(self.y[1:-1], corner.y) + 1)
        ind2 = (bisect.bisect_right(self.y[1:-1], corner.yy) + 1) + 1
        self.y.insert(ind1, corner.y)
        self.y.insert(ind2, corner.yy)
        self.n.insert(ind1, self.n[ind1-1])
        self.n.insert(ind2-1, self.n[ind2-1])
        for ind in range(ind1,ind2):
            self.n[ind] = self.n[ind] + 1

    def merge(self, corner): 
        ind1 = self.y[1:-1].index(corner.yy) + 1
        ind2 = self.y[1:-1].index(corner.y) + 1
        del self.y[ind1], self.y[ind2-1]
        for ind in range(ind1,ind2):
            self.n[ind] = self.n[ind] - 1
        del self.n[ind1], self.n[ind2-1]
        
def infcic(points, widths):

    npts, nscales = len(points), len(widths)
    Pn = np.empty((npts + 1, nscales))
            
    for k in range(nscales):
        
        width = widths[k]
        
        print('doing cell %s of %s ... (%.3f deg)' % (k + 1, nscales, width))
        
        # define inner boundary
        xmin, xmax = points[:,0].min() + 0.5 * width, points[:,0].max() - 0.5 * width
        ymin, ymax = points[:,1].min() + 0.5 * width, points[:,1].max() - 0.5 * width
        
        # place a cell (defined by two corners) over each point
        x1, y1 = (points[:,0] - 0.5 * width, points[:,1] - 0.5 * width) # lower left corner
        x2, y2 = (points[:,0] + 0.5 * width, points[:,1] + 0.5 * width) # upper right corner

        # reject cells that are fully outside inner boundary
        outside = ( (x1 >= xmax) | (y1 >= ymax) | (x2 <= xmin) | (y2 <= ymin) )
        x1, y1, x2, y2 = x1[~outside], y1[~outside], x2[~outside], y2[~outside]

        # clip cells that are partially outside inner boundary
        partial = ( (x1 < xmin) | (y1 < ymin) | (x2 > xmax) | (y2 > ymax) )
        x1[partial], y1[partial] = np.clip(x1[partial], xmin, xmax), np.clip(y1[partial], ymin, ymax)
        x2[partial], y2[partial] = np.clip(x2[partial], xmin, xmax), np.clip(y2[partial], ymin, ymax)

        # put cell corners into a list, sorted by ascending x-value
        ncells = len(x1)
        corners = []
        for j in range(ncells):
            corners.extend([Corner('L', x1[j], y1[j], y2[j]),Corner('R', x2[j], y2[j], y1[j])])
        corners = sorted(corners, key=lambda k: k.x)

        # sweep; update probabilities; if contact left (right) side of a cell, two markers added (dropped)
        sweep = Sweep(xmin, ymin, ymax)
        P = np.zeros(npts + 1)        
        for i in range(2*ncells):
            thiscorner = corners[i]
            if sweep.x != thiscorner.x:
                dx = thiscorner.x - sweep.x
                sweep.update(thiscorner.x)
                segments = np.diff(sweep.y)
                nmax = max(sweep.n)
                P[0:nmax+1] += np.bincount(sweep.n, segments * dx, minlength=nmax)
            if thiscorner.side == 'L':
                sweep.split(thiscorner)
            else:
                sweep.merge(thiscorner)        
                
        Pn[:,k] = P / P.sum()
        
    return Pn
