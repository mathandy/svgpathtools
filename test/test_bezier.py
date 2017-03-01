from __future__ import division, absolute_import, print_function
import numpy as np
import unittest
from svgpathtools.bezier import *
from svgpathtools.path import bpoints2bezier


class HigherOrderBezier:
    def __init__(self, bpoints):
        self.bpts = bpoints

    def bpoints(self):
        return self.bpts

    def point(self, t):
        return bezier_point(self.bpoints(), t)

    def __repr__(self):
        return str(self.bpts)


def random_polynomial(degree):
    return np.poly1d(np.random.rand(degree + 1))


def random_bezier(degree):
    if degree <= 3:
        return bpoints2bezier(polynomial2bezier(np.random.rand(degree + 1)))
    else:
        return HigherOrderBezier(np.random.rand(degree + 1))


class TestBezier2Polynomial(unittest.TestCase):
    def test_bezier2polynomial(self):
        tvals = np.linspace(0, 1, 10)
        for d in range(1, 10):
            b = random_bezier(d)
            p = np.poly1d(bezier2polynomial(b.bpoints()))
            for t in tvals:
                msg = ("degree {}\nt = {}\nb(t) = {}\n = {}\np(t) = \n{}\n = {}"
                       "".format(d, t, b, b.point(t), p, p(t)))
                self.assertAlmostEqual(b.point(t), p(t), msg=msg)


class TestPolynomial2Bezier(unittest.TestCase):
    def test_polynomial2bezier(self):
        tvals = np.linspace(0, 1, 10)
        for d in range(1, 3):
            p = random_polynomial(d)
            b = HigherOrderBezier(polynomial2bezier(p))
            for t in tvals:
                msg = ("degree {}\nt = {}\nb(t) = {}\n = {}\np(t) = \n{}\n = {}"
                       "".format(d, t, b, b.point(t), p, p(t)))
                self.assertAlmostEqual(b.point(t), p(t), msg=msg)


if __name__ == '__main__':
    unittest.main()
