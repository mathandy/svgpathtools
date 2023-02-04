from __future__ import division, absolute_import, print_function
import unittest
from svgpathtools import Document
from io import StringIO
from io import open  # overrides build-in open for compatibility with python2
from os.path import join, dirname
from sys import version_info


class TestDocument(unittest.TestCase):
    def test_from_file_path_string(self):
        """Test reading svg from file provided as path"""
        doc = Document(join(dirname(__file__), 'polygons.svg'))

        self.assertEqual(len(doc.paths()), 2)

    def test_from_file_path(self):
        """Test reading svg from file provided as path"""
        if version_info >= (3, 6):
            import pathlib
            doc = Document(pathlib.Path(__file__).parent / 'polygons.svg')

            self.assertEqual(len(doc.paths()), 2)

    def test_from_file_object(self):
        """Test reading svg from file object that has already been opened"""
        with open(join(dirname(__file__), 'polygons.svg'), 'r') as file:
            doc = Document(file)

            self.assertEqual(len(doc.paths()), 2)

    def test_from_stringio(self):
        """Test reading svg object contained in a StringIO object"""
        with open(join(dirname(__file__), 'polygons.svg'),
                  'r', encoding='utf-8') as file:
            # read entire file into string
            file_content = file.read()
            # prepare stringio object
            file_as_stringio = StringIO(file_content)

            doc = Document(file_as_stringio)

            self.assertEqual(len(doc.paths()), 2)

    def test_from_string(self):
        """Test reading svg object contained in a string"""
        with open(join(dirname(__file__), 'polygons.svg'),
                  'r', encoding='utf-8') as file:
            # read entire file into string
            file_content = file.read()

            doc = Document.from_svg_string(file_content)

            self.assertEqual(len(doc.paths()), 2)
