from __future__ import division, absolute_import, print_function
import unittest
from svgpathtools import Path, Line, Arc, svg2paths, svgstr2paths
from io import StringIO
from io import open  # overrides build-in open for compatibility with python2
import os
from os.path import join, dirname
from sys import version_info
import tempfile
import shutil

from svgpathtools.svg_to_paths import rect2pathd


class TestSVG2Paths(unittest.TestCase):
    def test_svg2paths_polygons(self):

        paths, _ = svg2paths(join(dirname(__file__), 'polygons.svg'))

        # triangular polygon test
        path = paths[0]
        path_correct = Path(Line(55.5+0j, 55.5+50j), 
                            Line(55.5+50j, 105.5+50j), 
                            Line(105.5+50j, 55.5+0j)
                            )
        self.assertTrue(path.isclosed())
        self.assertTrue(len(path)==3)
        self.assertTrue(path==path_correct)

        # triangular quadrilateral (with a redundant 4th "closure" point)
        path = paths[1]
        path_correct = Path(Line(0+0j, 0-100j),
                            Line(0-100j, 0.1-100j),
                            Line(0.1-100j, 0+0j),
                            Line(0+0j, 0+0j)  # result of redundant point
                            )
        self.assertTrue(path.isclosed())
        self.assertTrue(len(path)==4)
        self.assertTrue(path==path_correct)

    def test_svg2paths_ellipses(self):

        paths, _ = svg2paths(join(dirname(__file__), 'ellipse.svg'))

        # ellipse tests
        path_ellipse = paths[0]
        path_ellipse_correct = Path(Arc(50+100j, 50+50j, 0.0, True, False, 150+100j),
                                    Arc(150+100j, 50+50j, 0.0, True, False, 50+100j))
        self.assertTrue(len(path_ellipse)==2)
        self.assertTrue(path_ellipse==path_ellipse_correct)
        self.assertTrue(path_ellipse.isclosed())

        # circle tests
        paths, _ = svg2paths(join(dirname(__file__), 'circle.svg'))

        path_circle = paths[0]
        path_circle_correct = Path(Arc(50+100j, 50+50j, 0.0, True, False, 150+100j),
                                    Arc(150+100j, 50+50j, 0.0, True, False, 50+100j))
        self.assertTrue(len(path_circle)==2)
        self.assertTrue(path_circle==path_circle_correct)
        self.assertTrue(path_circle.isclosed())

        # test for issue #198 (circles not being closed)
        svg = u"""<?xml version="1.0" encoding="UTF-8"?>
        <svg xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" width="40mm" height="40mm" 
          viewBox="0 0 40 40" version="1.1">
            
          <g id="layer">
            <circle id="c1" cx="20.000" cy="20.000" r="11.000" />
            <circle id="c2" cx="20.000" cy="20.000" r="5.15" />
          </g>
        </svg>"""
        tmpdir = tempfile.mkdtemp()
        svgfile = os.path.join(tmpdir, 'test.svg')
        with open(svgfile, 'w') as f:
            f.write(svg)
        paths, _ = svg2paths(svgfile)
        self.assertEqual(len(paths), 2)
        self.assertTrue(paths[0].isclosed())
        self.assertTrue(paths[1].isclosed())
        shutil.rmtree(tmpdir)

    def test_rect2pathd(self):
        non_rounded_dict = {"x": "10", "y": "10", "width": "100", "height": "100"}
        self.assertEqual(
            rect2pathd(non_rounded_dict),
            "M10.0 10.0 L 110.0 10.0 L 110.0 110.0 L 10.0 110.0 z",
        )

        non_rounded_svg = """<?xml version="1.0" encoding="UTF-8"?>
          <svg xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" width="200mm" height="200mm" version="1.1">
              <rect id="non_rounded" x="10" y="10" width="100" height="100" />
          </svg>"""

        paths, _ = svg2paths(StringIO(non_rounded_svg))
        self.assertEqual(len(paths), 1)
        self.assertTrue(paths[0].isclosed())
        self.assertEqual(
            paths[0].d(use_closed_attrib=True),
            "M 10.0,10.0 L 110.0,10.0 L 110.0,110.0 L 10.0,110.0 Z",
        )
        self.assertEqual(
            paths[0].d(use_closed_attrib=False),
            "M 10.0,10.0 L 110.0,10.0 L 110.0,110.0 L 10.0,110.0 L 10.0,10.0",
        )

        rounded_dict = {"x": "10", "y": "10", "width": "100","height": "100", "rx": "15", "ry": "12"}
        self.assertEqual(
            rect2pathd(rounded_dict),
            "M 25.0 10.0 L 95.0 10.0 A 15.0 12.0 0 0 1 110.0 22.0 L 110.0 98.0 A 15.0 12.0 0 0 1 95.0 110.0 L 25.0 110.0 A 15.0 12.0 0 0 1 10.0 98.0 L 10.0 22.0 A 15.0 12.0 0 0 1 25.0 10.0 z",
        )

        rounded_svg = """<?xml version="1.0" encoding="UTF-8"?>
          <svg xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" width="200mm" height="200mm" version="1.1">
              <rect id="rounded" x="10" y="10" width="100" height ="100" rx="15" ry="12" />
          </svg>"""

        paths, _ = svg2paths(StringIO(rounded_svg))
        self.assertEqual(len(paths), 1)
        self.assertTrue(paths[0].isclosed())
        self.assertEqual(
            paths[0].d(),
            "M 25.0,10.0 L 95.0,10.0 A 15.0,12.0 0.0 0,1 110.0,22.0 L 110.0,98.0 A 15.0,12.0 0.0 0,1 95.0,110.0 L 25.0,110.0 A 15.0,12.0 0.0 0,1 10.0,98.0 L 10.0,22.0 A 15.0,12.0 0.0 0,1 25.0,10.0",
        )

    def test_from_file_path_string(self):
        """Test reading svg from file provided as path"""
        paths, _ = svg2paths(join(dirname(__file__), 'polygons.svg'))

        self.assertEqual(len(paths), 2)

    def test_from_file_path(self):
        """Test reading svg from file provided as pathlib POSIXPath"""
        if version_info >= (3, 6):
            import pathlib
            paths, _ = svg2paths(pathlib.Path(__file__).parent / 'polygons.svg')

            self.assertEqual(len(paths), 2)

    def test_from_file_object(self):
        """Test reading svg from file object that has already been opened"""
        with open(join(dirname(__file__), 'polygons.svg'), 'r') as file:
            paths, _ = svg2paths(file)

            self.assertEqual(len(paths), 2)

    def test_from_stringio(self):
        """Test reading svg object contained in a StringIO object"""
        with open(join(dirname(__file__), 'polygons.svg'),
                  'r', encoding='utf-8') as file:
            # read entire file into string
            file_content = file.read()
            # prepare stringio object
            file_as_stringio = StringIO(file_content)

            paths, _ = svg2paths(file_as_stringio)

            self.assertEqual(len(paths), 2)

    def test_from_string(self):
        """Test reading svg object contained in a string"""
        with open(join(dirname(__file__), 'polygons.svg'),
                  'r', encoding='utf-8') as file:
            # read entire file into string
            file_content = file.read()

            paths, _ = svgstr2paths(file_content)

            self.assertEqual(len(paths), 2)

    def test_svg2paths_polygon_no_points(self):

        paths, _ = svg2paths(join(dirname(__file__), 'polygons_no_points.svg'))

        path = paths[0]
        path_correct = Path()
        self.assertTrue(len(path)==0)
        self.assertTrue(path==path_correct)

        path = paths[1]
        self.assertTrue(len(path)==0)
        self.assertTrue(path==path_correct)

    def test_svg2paths_polyline_tests(self):

        paths, _ = svg2paths(join(dirname(__file__), 'polyline.svg'))

        path = paths[0]
        path_correct = Path(Line(59+185j, 98+203j),
                            Line(98+203j, 108+245j),
                            Line(108+245j, 82+279j),
                            Line(82+279j, 39+280j),
                            Line(39+280j, 11+247j),
                            Line(11+247j, 19+205j))
        self.assertFalse(path.isclosed())
        self.assertTrue(len(path)==6)
        self.assertTrue(path==path_correct)

        path = paths[1]
        path_correct = Path(Line(220+50j, 267+84j),
                            Line(267+84j, 249+140j),
                            Line(249+140j, 190+140j),
                            Line(190+140j, 172+84j),
                            Line(172+84j, 220+50j))
        self.assertTrue(path.isclosed())
        self.assertTrue(len(path)==5)
        self.assertTrue(path==path_correct)


if __name__ == '__main__':
    unittest.main()
