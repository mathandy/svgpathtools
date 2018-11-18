from __future__ import division, absolute_import, print_function
import unittest
from svgpathtools import *

class TestPathD(unittest.TestCase):

    def test_d(self):
        # the following two path represent the same path but in absolute and relative forms
        abs_s = 'M 38.0,130.0 C 37.0,132.0 38.0,136.0 40.0,137.0 L 85.0,161.0 C 87.0,162.0 91.0,162.0 93.0,160.0 L 127.0,133.0 C 129.0,131.0 129.0,128.0 127.0,126.0 L 80.0,70.0 C 78.0,67.0 75.0,68.0 74.0,70.0 Z'
        rel_s = 'm 38.0,130.0 c -1.0,2.0 0.0,6.0 2.0,7.0 l 45.0,24.0 c 2.0,1.0 6.0,1.0 8.0,-1.0 l 34.0,-27.0 c 2.0,-2.0 2.0,-5.0 0.0,-7.0 l -47.0,-56.0 c -2.0,-3.0 -5.0,-2.0 -6.0,0.0 z'
        path1 = parse_path(abs_s)
        path2 = parse_path(rel_s)
        self.assertEqual(path1.d(use_closed_attrib=True), abs_s)
        self.assertEqual(path2.d(use_closed_attrib=True), abs_s)
        self.assertEqual(path1.d(use_closed_attrib=True, rel=True), rel_s)
        self.assertEqual(path2.d(use_closed_attrib=True, rel=True), rel_s)

if __name__ == '__main__':
    unittest.main()

