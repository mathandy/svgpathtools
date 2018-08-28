from __future__ import division, absolute_import, print_function
import unittest
from svgpathtools import *
from os.path import join, dirname


class TestSaxGroups(unittest.TestCase):

    def check_values(self, v, z):
        # Check that the components of 2D vector v match the components
        # of complex number z
        self.assertAlmostEqual(v[0], z.real)
        self.assertAlmostEqual(v[1], z.imag)

    def test_parse_display(self):
        doc = SaxDocument(join(dirname(__file__), 'groups.svg'))
        doc.display()
