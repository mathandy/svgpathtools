"""The goal of this gist is to show how to compute many points on a path 
quickly using NumPy arrays.  I.e. there's a much faster way than using, say
[some_path.point(t) for t in many_tvals].  The example below assumes the 
`Path` object is composed entirely of `CubicBezier` objects, but this can
easily be generalized to paths containing `Line` and `QuadraticBezier` objects
also.  
Note: The relevant matrix transformation for quadratics can be found in the
svgpathtools.bezier module."""
from __future__ import print_function
import numpy as np
from svgpathtools import *


class HigherOrderBezier:
    def __init__(self, bpoints):
        self.bpts = bpoints

    def bpoints(self):
        return self.bpts

    def point(self, t):
        return bezier_point(self.bpoints(), t)

    def __repr__(self):
        return str(self.bpts)


def random_bezier(degree):
    if degree <= 3:
        return bpoints2bezier(polynomial2bezier(np.random.rand(degree + 1)))
    else:
        return HigherOrderBezier(np.random.rand(degree + 1))


def points_in_each_seg_slow(path, tvals):
    return [seg.poly()(tvals) for seg in path]


def points_in_each_seg(path, tvals):
    """Compute seg.point(t) for each seg in path and each t in tvals."""
    A = np.array([[-1,  3, -3,  1], # transforms cubic bez to standard poly
                  [ 3, -6,  3,  0],
                  [-3,  3,  0,  0],
                  [ 1,  0,  0,  0]])
    B = [seg.bpoints() for seg in path]
    return np.dot(B, np.dot(A, np.power(tvals, [[3],[2],[1],[0]])))


if __name__ == '__main__':
    num_segs = 1000
    testpath = Path(*[random_bezier(3) for dummy in range(num_segs)])
    tvals = np.linspace(0, 1, 10)

    pts = points_in_each_seg(testpath, tvals)
    pts_check = points_in_each_seg_slow(testpath, tvals)
    print(np.max(pts - pts_check))
