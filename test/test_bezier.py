from __future__ import division, absolute_import, print_function
import numpy as np
import unittest
from svgpathtools.bezier import bezier_point, bezier2polynomial, polynomial2bezier, bezier_bounding_box, bezier_real_minmax
from svgpathtools.path import bpoints2bezier, CubicBezier


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


class TestBezierBoundingBox(unittest.TestCase):
    def test_bezier_bounding_box(self):
        # This bezier curve has denominator == 0 but due to floating point arithmetic error it is not exactly 0
        zero_denominator_bezier_curve = CubicBezier(612.547 + 109.3261j, 579.967 - 19.4422j, 428.0344 - 19.4422j, 395.4374 + 109.3261j)
        zero_denom_xmin, zero_denom_xmax, zero_denom_ymin, zero_denom_ymax = bezier_bounding_box(zero_denominator_bezier_curve)
        self.assertAlmostEqual(zero_denom_xmin, 395.437400, 5)
        self.assertAlmostEqual(zero_denom_xmax, 612.547, 5)
        self.assertAlmostEqual(zero_denom_ymin, 12.7498749, 5)
        self.assertAlmostEqual(zero_denom_ymax, 109.3261, 5)

        # This bezier curve has global extrema at the start and end points
        start_end_bbox_bezier_curve = CubicBezier(886.8238 + 354.8439j, 884.4765 + 340.5983j, 877.6258 + 330.0518j, 868.2909 + 323.2453j)
        start_end_xmin, start_end_xmax, start_end_ymin, start_end_ymax = bezier_bounding_box(start_end_bbox_bezier_curve)
        self.assertAlmostEqual(start_end_xmin, 868.2909, 5)
        self.assertAlmostEqual(start_end_xmax, 886.8238, 5)
        self.assertAlmostEqual(start_end_ymin, 323.2453, 5)
        self.assertAlmostEqual(start_end_ymax, 354.8439, 5)

        # This bezier curve is to cover some random case where at least one of the global extrema is not the start or end point
        general_bezier_curve = CubicBezier(295.2282 + 402.0233j, 310.3734 + 355.5329j, 343.547 + 340.5983j, 390.122 + 355.7018j)
        general_xmin, general_xmax, general_ymin, general_ymax = bezier_bounding_box(general_bezier_curve)
        self.assertAlmostEqual(general_xmin, 295.2282, 5)
        self.assertAlmostEqual(general_xmax, 390.121999999, 5)
        self.assertAlmostEqual(general_ymin, 350.030030142, 5)
        self.assertAlmostEqual(general_ymax, 402.0233, 5)


class TestBezierRealMinMax(unittest.TestCase):
    def test_bezier_real_minmax(self):
        # This bezier curve has denominator == 0 but due to floating point arithmetic error it is not exactly 0
        zero_denominator_bezier_curve = [109.3261, -19.4422, -19.4422, 109.3261]
        zero_denominator_minmax = bezier_real_minmax(zero_denominator_bezier_curve)
        self.assertAlmostEqual(zero_denominator_minmax[0], 12.7498749, 5)
        self.assertAlmostEqual(zero_denominator_minmax[1], 109.3261, 5)

        # This bezier curve has global extrema at the start and end points
        start_end_bbox_bezier_curve = [354.8439, 340.5983, 330.0518, 323.2453]
        start_end_bbox_minmax = bezier_real_minmax(start_end_bbox_bezier_curve)
        self.assertAlmostEqual(start_end_bbox_minmax[0], 323.2453, 5)
        self.assertAlmostEqual(start_end_bbox_minmax[1], 354.8439, 5)

        # This bezier curve is to cover some random case where at least one of the global extrema is not the start or end point
        general_bezier_curve = [402.0233, 355.5329, 340.5983, 355.7018]
        general_minmax = bezier_real_minmax(general_bezier_curve)
        self.assertAlmostEqual(general_minmax[0], 350.030030142, 5)
        self.assertAlmostEqual(general_minmax[1], 402.0233, 5)


if __name__ == '__main__':
    unittest.main()
