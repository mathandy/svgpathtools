# External dependencies
from __future__ import division, absolute_import, print_function
import unittest
import numpy as np

# Internal dependencies
from svgpathtools import *


class Test_polytools(unittest.TestCase):
    # def test_poly_roots(self):
    #     self.fail()

    def test_rational_limit(self):

        # (3x^3 + x)/(4x^2 - 2x) -> -1/2 as x->0
        f = np.poly1d([3, 0, 1, 0])
        g = np.poly1d([4, -2, 0])
        lim = rational_limit(f, g, 0)
        self.assertAlmostEqual(lim, -0.5)

        # (3x^2)/(4x^2 - 2x) -> 0 as x->0
        f = np.poly1d([3, 0, 0])
        g = np.poly1d([4, -2, 0])
        lim = rational_limit(f, g, 0)
        self.assertAlmostEqual(lim, 0)


if __name__ == '__main__':
    unittest.main()
