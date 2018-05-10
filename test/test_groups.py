from __future__ import division, absolute_import, print_function
import unittest
from svgpathtools import *
from os.path import join, dirname


class TestGroups(unittest.TestCase):

    def test_group_flatten(self):
        doc = Document(join(dirname(__file__), 'groups.svg'))

        result = doc.flatten_all_paths()
        print('\nnumber of paths: '+str(len(result)))
        self.assertGreater(len(result), 0)
