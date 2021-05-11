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
        doc = SaxDocument(join(dirname(__file__), 'transforms.svg'))
        # doc.display()
        for i, node in enumerate(doc.tree):
            values = node
            path_value = values['d']
            matrix = values['matrix']
            self.assertTrue(values is not None)
            self.assertTrue(path_value is not None)
            if i == 0:
                self.assertEqual(values['fill'], 'red')
            if i == 8 or i == 7:
                self.assertEqual(matrix, None)
            if i == 9:
                self.assertEqual(values['fill'], 'lime')
