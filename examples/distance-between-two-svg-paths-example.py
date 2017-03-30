from svgpathtools import *

# create some example paths
path1 = CubicBezier(1,2+3j,3-5j,4+1j)
path2 = path1.rotated(60).translated(3)

# find minimizer
from scipy.optimize import fminbound
def dist(t):
	return path1.radialrange(path2.point(t))[0][0]
T2 = fminbound(dist, 0, 1)

# Let's do a visual check
pt2 = path2.point(T2)
T1 = path1.radialrange(pt2)[0][1]
pt1 = path1.point(T1)
disvg([path1, path2, Line(pt1, pt2)], 'grb', nodes=[pt1, pt2])
