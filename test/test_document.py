from __future__ import division, absolute_import, print_function
import unittest
from svgpathtools import *
from io import StringIO
from os.path import join, dirname

class TestDocument(unittest.TestCase):
    def test_from_file_path(self):
        """ Test reading svg from file provided as path """
        doc = Document(join(dirname(__file__), 'polygons.svg'))

        self.assertEqual(len(doc.paths()), 2)

    def test_from_file_object(self):
        """ Test reading svg from file object that has already been opened """
        with open(join(dirname(__file__), 'polygons.svg'), 'r') as file:
            doc = Document(file)

            self.assertEqual(len(doc.paths()), 2)

    def test_from_stringio(self):
        """ Test reading svg object contained in a StringIO object """
        with open(join(dirname(__file__), 'polygons.svg'), 'r') as file:
            # read entire file into string
            file_content = file.read()
            # prepare stringio object
            file_as_stringio = StringIO()
            # paste file content into it
            file_as_stringio.write(file_content)
            # reset curser to its beginning
            file_as_stringio.seek(0)

            doc = Document(file_as_stringio)

            self.assertEqual(len(doc.paths()), 2)

    def test_from_string_without_svg_attrs(self):
        """ Test reading svg object contained in a string without svg attributes"""
        with open(join(dirname(__file__), 'polygons.svg'), 'r') as file:
            # read entire file into string
            file_content = file.read()

            doc = Document.from_svg_string(file_content)

            self.assertEqual(len(doc.paths()), 2)
