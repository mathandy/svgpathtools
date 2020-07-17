# Note: This file was taken mostly as is from the svg.path module (v 2.0)
from __future__ import division, absolute_import, print_function
import unittest
from svgpathtools import *
import svgpathtools
import numpy as np


def construct_rotation_tf(a, x, y):
    a = a * np.pi / 180.0
    tf_offset = np.identity(3)
    tf_offset[0:2, 2:3] = np.array([[x], [y]])
    tf_rotate = np.identity(3)
    tf_rotate[0:2, 0:2] = np.array([[np.cos(a), -np.sin(a)],
                                    [np.sin(a), np.cos(a)]])
    tf_offset_neg = np.identity(3)
    tf_offset_neg[0:2, 2:3] = np.array([[-x], [-y]])

    return tf_offset.dot(tf_rotate).dot(tf_offset_neg)


class TestParser(unittest.TestCase):

    def test_svg_examples(self):
        """Examples from the SVG spec"""
        path1 = parse_path('M 100 100 L 300 100 L 200 300 z')
        self.assertEqual(path1, Path(Line(100 + 100j, 300 + 100j),
                                     Line(300 + 100j, 200 + 300j),
                                     Line(200 + 300j, 100 + 100j)))
        self.assertTrue(path1.isclosed())

        # for Z command behavior when there is multiple subpaths
        path1 = parse_path('M 0 0 L 50 20 M 100 100 L 300 100 L 200 300 z')
        self.assertEqual(path1, Path(Line(0 + 0j, 50 + 20j),
                                     Line(100 + 100j, 300 + 100j),
                                     Line(300 + 100j, 200 + 300j),
                                     Line(200 + 300j, 100 + 100j)))

        path1 = parse_path('M 100 100 L 200 200')
        path2 = parse_path('M100 100L200 200')
        self.assertEqual(path1, path2)

        path1 = parse_path('M 100 200 L 200 100 L -100 -200')
        path2 = parse_path('M 100 200 L 200 100 -100 -200')
        self.assertEqual(path1, path2)

        path1 = parse_path("""M100,200 C100,100 250,100 250,200
                              S400,300 400,200""")
        self.assertEqual(path1, Path(CubicBezier(100 + 200j,
                                                 100 + 100j,
                                                 250 + 100j,
                                                 250 + 200j),
                                     CubicBezier(250 + 200j,
                                                 250 + 300j,
                                                 400 + 300j,
                                                 400 + 200j)))

        path1 = parse_path('M100,200 C100,100 400,100 400,200')
        self.assertEqual(path1, Path(CubicBezier(100 + 200j,
                                                 100 + 100j,
                                                 400 + 100j,
                                                 400 + 200j)))

        path1 = parse_path('M100,500 C25,400 475,400 400,500')
        self.assertEqual(path1, Path(CubicBezier(100 + 500j,
                                                 25 + 400j,
                                                 475 + 400j,
                                                 400 + 500j)))

        path1 = parse_path('M100,800 C175,700 325,700 400,800')
        self.assertEqual(path1, Path(CubicBezier(100 + 800j,
                                                 175 + 700j,
                                                 325 + 700j,
                                                 400 + 800j)))

        path1 = parse_path('M600,200 C675,100 975,100 900,200')
        self.assertEqual(path1, Path(CubicBezier(600 + 200j,
                                                 675 + 100j,
                                                 975 + 100j,
                                                 900 + 200j)))

        path1 = parse_path('M600,500 C600,350 900,650 900,500')
        self.assertEqual(path1, Path(CubicBezier(600 + 500j,
                                                 600 + 350j,
                                                 900 + 650j,
                                                 900 + 500j)))

        path1 = parse_path("""M600,800 C625,700 725,700 750,800
                              S875,900 900,800""")
        self.assertEqual(path1, Path(CubicBezier(600 + 800j,
                                                 625 + 700j,
                                                 725 + 700j,
                                                 750 + 800j),
                                     CubicBezier(750 + 800j,
                                                 775 + 900j,
                                                 875 + 900j,
                                                 900 + 800j)))

        path1 = parse_path('M200,300 Q400,50 600,300 T1000,300')
        self.assertEqual(path1, Path(QuadraticBezier(200 + 300j,
                                                     400 + 50j,
                                                     600 + 300j),
                                     QuadraticBezier(600 + 300j,
                                                     800 + 550j,
                                                     1000 + 300j)))

        path1 = parse_path('M300,200 h-150 a150,150 0 1,0 150,-150 z')
        self.assertEqual(path1, Path(Line(300 + 200j, 150 + 200j),
                                Arc(150 + 200j, 150 + 150j, 0, 1, 0, 300 + 50j),
                                Line(300 + 50j, 300 + 200j)))

        path1 = parse_path('M275,175 v-150 a150,150 0 0,0 -150,150 z')
        self.assertEqual(path1,
                         Path(Line(275 + 175j, 275 + 25j),
                              Arc(275 + 25j, 150 + 150j, 0, 0, 0, 125 + 175j),
                              Line(125 + 175j, 275 + 175j)))

        path1 = parse_path("""M600,350 l 50,-25
                              a25,25 -30 0,1 50,-25 l 50,-25
                              a25,50 -30 0,1 50,-25 l 50,-25
                              a25,75 -30 0,1 50,-25 l 50,-25
                              a25,100 -30 0,1 50,-25 l 50,-25""")
        self.assertEqual(path1,
                         Path(Line(600 + 350j, 650 + 325j),
                              Arc(650 + 325j, 25 + 25j, -30, 0, 1, 700 + 300j),
                              Line(700 + 300j, 750 + 275j),
                              Arc(750 + 275j, 25 + 50j, -30, 0, 1, 800 + 250j),
                              Line(800 + 250j, 850 + 225j),
                              Arc(850 + 225j, 25 + 75j, -30, 0, 1, 900 + 200j),
                              Line(900 + 200j, 950 + 175j),
                              Arc(950 + 175j, 25 + 100j, -30, 0, 1, 1000 + 150j),
                              Line(1000 + 150j, 1050 + 125j)))

    def test_others(self):
        # Other paths that need testing:

        # Relative moveto:
        path1 = parse_path('M 0 0 L 50 20 m 50 80 L 300 100 L 200 300 z')
        self.assertEqual(path1, Path(Line(0 + 0j, 50 + 20j),
                                     Line(100 + 100j, 300 + 100j),
                                     Line(300 + 100j, 200 + 300j),
                                     Line(200 + 300j, 100 + 100j)))

        # Initial smooth and relative CubicBezier
        path1 = parse_path("""M100,200 s 150,-100 150,0""")
        self.assertEqual(path1,
                         Path(CubicBezier(100 + 200j,
                                          100 + 200j,
                                          250 + 100j,
                                          250 + 200j)))

        # Initial smooth and relative QuadraticBezier
        path1 = parse_path("""M100,200 t 150,0""")
        self.assertEqual(path1,
                         Path(QuadraticBezier(100 + 200j,
                                              100 + 200j,
                                              250 + 200j)))

        # Relative QuadraticBezier
        path1 = parse_path("""M100,200 q 0,0 150,0""")
        self.assertEqual(path1,
                         Path(QuadraticBezier(100 + 200j,
                                              100 + 200j,
                                              250 + 200j)))

    def test_negative(self):
        """You don't need spaces before a minus-sign"""
        path1 = parse_path('M100,200c10-5,20-10,30-20')
        path2 = parse_path('M 100 200 c 10 -5 20 -10 30 -20')
        self.assertEqual(path1, path2)

    def test_numbers(self):
        """Exponents and other number format cases"""
        # It can be e or E, the plus is optional, and a minimum of
        # +/-3.4e38 must be supported.
        path1 = parse_path('M-3.4e38 3.4E+38L-3.4E-38,3.4e-38')
        path2 = Path(Line(-3.4e+38 + 3.4e+38j, -3.4e-38 + 3.4e-38j))
        self.assertEqual(path1, path2)

    def test_errors(self):
        self.assertRaises(ValueError, parse_path,
                          'M 100 100 L 200 200 Z 100 200')


    def test_transform(self):

        tf_matrix = svgpathtools.parser.parse_transform(
            'matrix(1.0 2.0 3.0 4.0 5.0 6.0)')
        expected_tf_matrix = np.identity(3)
        expected_tf_matrix[0:2, 0:3] = np.array([[1.0, 3.0, 5.0],
                                                 [2.0, 4.0, 6.0]])
        self.assertTrue(np.array_equal(expected_tf_matrix, tf_matrix))

        # Try a test with no y specified
        expected_tf_translate = np.identity(3)
        expected_tf_translate[0, 2] = -36
        self.assertTrue(np.array_equal(
            expected_tf_translate,
            svgpathtools.parser.parse_transform('translate(-36)')
        ))

        # Now specify y
        expected_tf_translate[1, 2] = 45.5
        tf_translate = svgpathtools.parser.parse_transform(
            'translate(-36 45.5)')
        self.assertTrue(np.array_equal(expected_tf_translate, tf_translate))

        # Try a test with no y specified
        expected_tf_scale = np.identity(3)
        expected_tf_scale[0, 0] = 10
        expected_tf_scale[1, 1] = 10
        self.assertTrue(np.array_equal(
            expected_tf_scale,
            svgpathtools.parser.parse_transform('scale(10)')
        ))

        # Now specify y
        expected_tf_scale[1, 1] = 0.5
        tf_scale = svgpathtools.parser.parse_transform('scale(10 0.5)')
        self.assertTrue(np.array_equal(expected_tf_scale, tf_scale))

        tf_rotation = svgpathtools.parser.parse_transform('rotate(-10 50 100)')
        expected_tf_rotation = construct_rotation_tf(-10, 50, 100)
        self.assertTrue(np.array_equal(expected_tf_rotation, tf_rotation))

        # Try a test with no offset specified
        self.assertTrue(np.array_equal(
            construct_rotation_tf(50, 0, 0),
            svgpathtools.parser.parse_transform('rotate(50)')
        ))

        expected_tf_skewx = np.identity(3)
        expected_tf_skewx[0, 1] = np.tan(40.0 * np.pi/180.0)
        tf_skewx = svgpathtools.parser.parse_transform('skewX(40)')
        self.assertTrue(np.array_equal(expected_tf_skewx, tf_skewx))

        expected_tf_skewy = np.identity(3)
        expected_tf_skewy[1, 0] = np.tan(30.0 * np.pi / 180.0)
        tf_skewy = svgpathtools.parser.parse_transform('skewY(30)')
        self.assertTrue(np.array_equal(expected_tf_skewy, tf_skewy))

        self.assertTrue(np.array_equal(
            tf_rotation.dot(tf_translate).dot(tf_skewx).dot(tf_scale),
            svgpathtools.parser.parse_transform(
                """rotate(-10 50 100)
                   translate(-36 45.5)
                   skewX(40)
                   scale(10 0.5)""")
        ))

    def test_pathd_init(self):
        path0 = Path('')
        path1 = parse_path("M 100 100 L 300 100 L 200 300 z")
        path2 = Path("M 100 100 L 300 100 L 200 300 z")
        self.assertEqual(path1, path2)

        path1 = parse_path("m 100 100 L 300 100 L 200 300 z", current_pos=50+50j)
        path2 = Path("m 100 100 L 300 100 L 200 300 z")
        self.assertNotEqual(path1, path2)

        path1 = parse_path("m 100 100 L 300 100 L 200 300 z")
        path2 = Path("m 100 100 L 300 100 L 200 300 z", current_pos=50 + 50j)
        self.assertNotEqual(path1, path2)

        path1 = parse_path("m 100 100 L 300 100 L 200 300 z", current_pos=50 + 50j)
        path2 = Path("m 100 100 L 300 100 L 200 300 z", current_pos=50 + 50j)
        self.assertEqual(path1, path2)

        path1 = parse_path("m 100 100 L 300 100 L 200 300 z", 50+50j)
        path2 = Path("m 100 100 L 300 100 L 200 300 z")
        self.assertNotEqual(path1, path2)

        path1 = parse_path("m 100 100 L 300 100 L 200 300 z")
        path2 = Path("m 100 100 L 300 100 L 200 300 z", 50 + 50j)
        self.assertNotEqual(path1, path2)

        path1 = parse_path("m 100 100 L 300 100 L 200 300 z", 50 + 50j)
        path2 = Path("m 100 100 L 300 100 L 200 300 z", 50 + 50j)
        self.assertEqual(path1, path2)

    def test_issue_99(self):
        p = Path("M  100 250    S  200 200   200 250     300 300   300 250")
        self.assertEqual(p.d(useSandT=True), 'M 100.0,250.0 S 200.0,200.0 200.0,250.0 S 300.0,300.0 300.0,250.0')
        self.assertEqual(p.d(),
                         'M 100.0,250.0 C 100.0,250.0 200.0,200.0 200.0,250.0 C 200.0,300.0 300.0,300.0 300.0,250.0')
        self.assertNotEqual(p.d(),
                         'M 100.0,250.0 C 100.0,250.0 200.0,200.0 200.0,250.0 C 200.0,250.0 300.0,300.0 300.0,250.0')
