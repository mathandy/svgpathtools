from __future__ import division, absolute_import, print_function
from unittest import TestCase
from svgpathtools import *
from os.path import join, dirname

class TestSvg2pathsGroups(TestCase):
    def test_svg2paths(self):
        paths, _ = svg2paths(join(dirname(__file__), 'groups.svg'))

        # the paths should form crosses after being transformed
        self.assertTrue((len(paths) % 2) == 0)

        for i in range(len(paths)//2):
            self.assertTrue(len(paths[i * 2].intersect(paths[i * 2 + 1])) > 0, 'Path '+str(i * 2)+' does not intersect path '+str(i * 2 + 1)+'!')