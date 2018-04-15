# External dependencies
from __future__ import division, absolute_import, print_function
import unittest
from math import sqrt, pi
from operator import itemgetter
from numpy import poly1d, isclose
import random

# Internal dependencies
from svgpathtools import *

# An important note for those doing any debugging:
# ------------------------------------------------
# Most of these test points are not calculated separately, as that would
# take too long and be too error prone. Instead the curves have been verified
# to be correct visually with the disvg() function.


def random_line():
    x = (random.random() - 0.5) * 2000
    y = (random.random() - 0.5) * 2000
    start = complex(x, y)

    x = (random.random() - 0.5) * 2000
    y = (random.random() - 0.5) * 2000
    end = complex(x, y)

    return Line(start, end)


def random_arc():
    x = (random.random() - 0.5) * 2000
    y = (random.random() - 0.5) * 2000
    start = complex(x, y)

    x = (random.random() - 0.5) * 2000
    y = (random.random() - 0.5) * 2000
    end = complex(x, y)

    x = (random.random() - 0.5) * 2000
    y = (random.random() - 0.5) * 2000
    radius = complex(x, y)

    large_arc = random.choice([True, False])
    sweep = random.choice([True, False])

    return Arc(start=start, radius=radius, rotation=0.0, large_arc=large_arc, sweep=sweep, end=end)


def assert_intersections(a_seg, b_seg, intersections, count):
    if count != None:
        assert(len(intersections) == count)
    for i in intersections:
        assert(i[0] >= 0.0)
        assert(i[0] <= 1.0)
        assert(i[1] >= 0.0)
        assert(i[1] <= 1.0)
        assert(isclose(a_seg.point(i[0]), b_seg.point(i[1])))


class LineTest(unittest.TestCase):

    def test_lines(self):
        # These points are calculated, and not just regression tests.

        line1 = Line(0j, 400 + 0j)
        self.assertAlmostEqual(line1.point(0), 0j)
        self.assertAlmostEqual(line1.point(0.3), (120 + 0j))
        self.assertAlmostEqual(line1.point(0.5), (200 + 0j))
        self.assertAlmostEqual(line1.point(0.9), (360 + 0j))
        self.assertAlmostEqual(line1.point(1), (400 + 0j))
        self.assertAlmostEqual(line1.length(), 400)

        line2 = Line(400 + 0j, 400 + 300j)
        self.assertAlmostEqual(line2.point(0), (400 + 0j))
        self.assertAlmostEqual(line2.point(0.3), (400 + 90j))
        self.assertAlmostEqual(line2.point(0.5), (400 + 150j))
        self.assertAlmostEqual(line2.point(0.9), (400 + 270j))
        self.assertAlmostEqual(line2.point(1), (400 + 300j))
        self.assertAlmostEqual(line2.length(), 300)

        line3 = Line(400 + 300j, 0j)
        self.assertAlmostEqual(line3.point(0), (400 + 300j))
        self.assertAlmostEqual(line3.point(0.3), (280 + 210j))
        self.assertAlmostEqual(line3.point(0.5), (200 + 150j))
        self.assertAlmostEqual(line3.point(0.9), (40 + 30j))
        self.assertAlmostEqual(line3.point(1), 0j)
        self.assertAlmostEqual(line3.length(), 500)

    def test_equality(self):
        # This is to test the __eq__ and __ne__ methods, so we can't use
        # assertEqual and assertNotEqual
        line = Line(0j, 400 + 0j)
        self.assertTrue(line == Line(0, 400))
        self.assertTrue(line != Line(100, 400))
        self.assertFalse(line == str(line))
        self.assertTrue(line != str(line))
        self.assertFalse(
            CubicBezier(600 + 500j, 600 + 350j, 900 + 650j, 900 + 500j) == line)

    def test_point_to_t(self):
        l = Line(start=(0+0j), end=(0+10j))
        assert(l.point_to_t(0+0j) == 0.0)
        assert(isclose(l.point_to_t(0+5j), 0.5))
        assert(l.point_to_t(0+10j) == 1.0)
        assert(l.point_to_t(1+0j) == None)
        assert(l.point_to_t(0-1j) == None)
        assert(l.point_to_t(0+11j) == None)

        l = Line(start=(0+0j), end=(10+10j))
        assert(l.point_to_t(0+0j) == 0.0)
        assert(isclose(l.point_to_t(5+5j), 0.5))
        assert(l.point_to_t(10+10j) == 1.0)
        assert(l.point_to_t(1+0j) == None)
        assert(l.point_to_t(0-1j) == None)
        assert(l.point_to_t(0+11j) == None)
        assert(l.point_to_t(10.001+10.001j) == None)
        assert(l.point_to_t(-0.001-0.001j) == None)

        l = Line(start=(0+0j), end=(10+0j))
        assert(l.point_to_t(0+0j) == 0.0)
        assert(isclose(l.point_to_t(5+0j), 0.5))
        assert(l.point_to_t(10+0j) == 1.0)
        assert(l.point_to_t(0+1j) == None)
        assert(l.point_to_t(0-1j) == None)
        assert(l.point_to_t(0+11j) == None)
        assert(l.point_to_t(10.001+0j) == None)
        assert(l.point_to_t(-0.001-0j) == None)

        l = Line(start=(-2-1j), end=(11-20j))
        assert(l.point_to_t(-2-1j) == 0.0)
        assert(isclose(l.point_to_t(4.5-10.5j), 0.5))
        assert(l.point_to_t(11-20j) == 1.0)
        assert(l.point_to_t(0+1j) == None)
        assert(l.point_to_t(0-1j) == None)
        assert(l.point_to_t(0+11j) == None)
        assert(l.point_to_t(10.001+0j) == None)
        assert(l.point_to_t(-0.001-0j) == None)

        l = Line(start=(40.234-32.613j), end=(12.7-32.613j))
        assert(l.point_to_t(40.234-32.613j) == 0.0)
        assert(isclose(l.point_to_t(33.3505-32.613j), 0.25))
        assert(isclose(l.point_to_t(26.467-32.613j), 0.50))
        assert(isclose(l.point_to_t(19.5835-32.613j), 0.75))
        assert(l.point_to_t(12.7-32.613j) == 1.0)
        assert(l.point_to_t(40.25-32.613j) == None)
        assert(l.point_to_t(12.65-32.613j) == None)
        assert(l.point_to_t(11-20j) == None)
        assert(l.point_to_t(0+1j) == None)
        assert(l.point_to_t(0-1j) == None)
        assert(l.point_to_t(0+11j) == None)
        assert(l.point_to_t(10.001+0j) == None)
        assert(l.point_to_t(-0.001-0j) == None)

        random.seed()
        for line_index in range(100):
            l = random_line()
            print(l)
            for t_index in range(100):
                orig_t = random.random()
                p = l.point(orig_t)
                computed_t = l.point_to_t(p)
                print("orig_t=%f, p=%s, computed_t=%f" % (orig_t, p, computed_t))
                assert(isclose(orig_t, computed_t))




class CubicBezierTest(unittest.TestCase):
    def test_approx_circle(self):
        """This is a approximate circle drawn in Inkscape"""

        cub1 = CubicBezier(
            complex(0, 0),
            complex(0, 109.66797),
            complex(-88.90345, 198.57142),
            complex(-198.57142, 198.57142)
        )

        self.assertAlmostEqual(cub1.point(0), 0j)
        self.assertAlmostEqual(cub1.point(0.1), (-2.59896457 + 32.20931647j))
        self.assertAlmostEqual(cub1.point(0.2), (-10.12330256 + 62.76392816j))
        self.assertAlmostEqual(cub1.point(0.3), (-22.16418039 + 91.25500149j))
        self.assertAlmostEqual(cub1.point(0.4), (-38.31276448 + 117.27370288j))
        self.assertAlmostEqual(cub1.point(0.5), (-58.16022125 + 140.41119875j))
        self.assertAlmostEqual(cub1.point(0.6), (-81.29771712 + 160.25865552j))
        self.assertAlmostEqual(cub1.point(0.7), (-107.31641851 + 176.40723961j))
        self.assertAlmostEqual(cub1.point(0.8), (-135.80749184 + 188.44811744j))
        self.assertAlmostEqual(cub1.point(0.9), (-166.36210353 + 195.97245543j))
        self.assertAlmostEqual(cub1.point(1), (-198.57142 + 198.57142j))

        cub2 = CubicBezier(
            complex(-198.57142, 198.57142),
            complex(-109.66797 - 198.57142, 0 + 198.57142),
            complex(-198.57143 - 198.57142, -88.90345 + 198.57142),
            complex(-198.57143 - 198.57142, 0),
        )

        self.assertAlmostEqual(cub2.point(0), (-198.57142 + 198.57142j))
        self.assertAlmostEqual(cub2.point(0.1), (-230.78073675 + 195.97245543j))
        self.assertAlmostEqual(cub2.point(0.2), (-261.3353492 + 188.44811744j))
        self.assertAlmostEqual(cub2.point(0.3), (-289.82642365 + 176.40723961j))
        self.assertAlmostEqual(cub2.point(0.4), (-315.8451264 + 160.25865552j))
        self.assertAlmostEqual(cub2.point(0.5), (-338.98262375 + 140.41119875j))
        self.assertAlmostEqual(cub2.point(0.6), (-358.830082 + 117.27370288j))
        self.assertAlmostEqual(cub2.point(0.7), (-374.97866745 + 91.25500149j))
        self.assertAlmostEqual(cub2.point(0.8), (-387.0195464 + 62.76392816j))
        self.assertAlmostEqual(cub2.point(0.9), (-394.54388515 + 32.20931647j))
        self.assertAlmostEqual(cub2.point(1), (-397.14285 + 0j))

        cub3 = CubicBezier(
            complex(-198.57143 - 198.57142, 0),
            complex(0 - 198.57143 - 198.57142, -109.66797),
            complex(88.90346 - 198.57143 - 198.57142, -198.57143),
            complex(-198.57142, -198.57143)
        )

        self.assertAlmostEqual(cub3.point(0), (-397.14285 + 0j))
        self.assertAlmostEqual(cub3.point(0.1), (-394.54388515 - 32.20931675j))
        self.assertAlmostEqual(cub3.point(0.2), (-387.0195464 - 62.7639292j))
        self.assertAlmostEqual(cub3.point(0.3), (-374.97866745 - 91.25500365j))
        self.assertAlmostEqual(cub3.point(0.4), (-358.830082 - 117.2737064j))
        self.assertAlmostEqual(cub3.point(0.5), (-338.98262375 - 140.41120375j))
        self.assertAlmostEqual(cub3.point(0.6), (-315.8451264 - 160.258662j))
        self.assertAlmostEqual(cub3.point(0.7), (-289.82642365 - 176.40724745j))
        self.assertAlmostEqual(cub3.point(0.8), (-261.3353492 - 188.4481264j))
        self.assertAlmostEqual(cub3.point(0.9), (-230.78073675 - 195.97246515j))
        self.assertAlmostEqual(cub3.point(1), (-198.57142 - 198.57143j))

        cub4 = CubicBezier(
            complex(-198.57142, -198.57143),
            complex(109.66797 - 198.57142, 0 - 198.57143),
            complex(0, 88.90346 - 198.57143),
            complex(0, 0),
        )

        self.assertAlmostEqual(cub4.point(0), (-198.57142 - 198.57143j))
        self.assertAlmostEqual(cub4.point(0.1), (-166.36210353 - 195.97246515j))
        self.assertAlmostEqual(cub4.point(0.2), (-135.80749184 - 188.4481264j))
        self.assertAlmostEqual(cub4.point(0.3), (-107.31641851 - 176.40724745j))
        self.assertAlmostEqual(cub4.point(0.4), (-81.29771712 - 160.258662j))
        self.assertAlmostEqual(cub4.point(0.5), (-58.16022125 - 140.41120375j))
        self.assertAlmostEqual(cub4.point(0.6), (-38.31276448 - 117.2737064j))
        self.assertAlmostEqual(cub4.point(0.7), (-22.16418039 - 91.25500365j))
        self.assertAlmostEqual(cub4.point(0.8), (-10.12330256 - 62.7639292j))
        self.assertAlmostEqual(cub4.point(0.9), (-2.59896457 - 32.20931675j))
        self.assertAlmostEqual(cub4.point(1), 0j)

    def test_svg_examples(self):

        # M100,200 C100,100 250,100 250,200
        path1 = CubicBezier(100 + 200j, 100 + 100j, 250 + 100j, 250 + 200j)
        self.assertAlmostEqual(path1.point(0), (100 + 200j))
        self.assertAlmostEqual(path1.point(0.3), (132.4 + 137j))
        self.assertAlmostEqual(path1.point(0.5), (175 + 125j))
        self.assertAlmostEqual(path1.point(0.9), (245.8 + 173j))
        self.assertAlmostEqual(path1.point(1), (250 + 200j))

        # S400,300 400,200
        path2 = CubicBezier(250 + 200j, 250 + 300j, 400 + 300j, 400 + 200j)
        self.assertAlmostEqual(path2.point(0), (250 + 200j))
        self.assertAlmostEqual(path2.point(0.3), (282.4 + 263j))
        self.assertAlmostEqual(path2.point(0.5), (325 + 275j))
        self.assertAlmostEqual(path2.point(0.9), (395.8 + 227j))
        self.assertAlmostEqual(path2.point(1), (400 + 200j))

        # M100,200 C100,100 400,100 400,200
        path3 = CubicBezier(100 + 200j, 100 + 100j, 400 + 100j, 400 + 200j)
        self.assertAlmostEqual(path3.point(0), (100 + 200j))
        self.assertAlmostEqual(path3.point(0.3), (164.8 + 137j))
        self.assertAlmostEqual(path3.point(0.5), (250 + 125j))
        self.assertAlmostEqual(path3.point(0.9), (391.6 + 173j))
        self.assertAlmostEqual(path3.point(1), (400 + 200j))

        # M100,500 C25,400 475,400 400,500
        path4 = CubicBezier(100 + 500j, 25 + 400j, 475 + 400j, 400 + 500j)
        self.assertAlmostEqual(path4.point(0), (100 + 500j))
        self.assertAlmostEqual(path4.point(0.3), (145.9 + 437j))
        self.assertAlmostEqual(path4.point(0.5), (250 + 425j))
        self.assertAlmostEqual(path4.point(0.9), (407.8 + 473j))
        self.assertAlmostEqual(path4.point(1), (400 + 500j))

        # M100,800 C175,700 325,700 400,800
        path5 = CubicBezier(100 + 800j, 175 + 700j, 325 + 700j, 400 + 800j)
        self.assertAlmostEqual(path5.point(0), (100 + 800j))
        self.assertAlmostEqual(path5.point(0.3), (183.7 + 737j))
        self.assertAlmostEqual(path5.point(0.5), (250 + 725j))
        self.assertAlmostEqual(path5.point(0.9), (375.4 + 773j))
        self.assertAlmostEqual(path5.point(1), (400 + 800j))

        # M600,200 C675,100 975,100 900,200
        path6 = CubicBezier(600 + 200j, 675 + 100j, 975 + 100j, 900 + 200j)
        self.assertAlmostEqual(path6.point(0), (600 + 200j))
        self.assertAlmostEqual(path6.point(0.3), (712.05 + 137j))
        self.assertAlmostEqual(path6.point(0.5), (806.25 + 125j))
        self.assertAlmostEqual(path6.point(0.9), (911.85 + 173j))
        self.assertAlmostEqual(path6.point(1), (900 + 200j))

        # M600,500 C600,350 900,650 900,500
        path7 = CubicBezier(600 + 500j, 600 + 350j, 900 + 650j, 900 + 500j)
        self.assertAlmostEqual(path7.point(0), (600 + 500j))
        self.assertAlmostEqual(path7.point(0.3), (664.8 + 462.2j))
        self.assertAlmostEqual(path7.point(0.5), (750 + 500j))
        self.assertAlmostEqual(path7.point(0.9), (891.6 + 532.4j))
        self.assertAlmostEqual(path7.point(1), (900 + 500j))

        # M600,800 C625,700 725,700 750,800
        path8 = CubicBezier(600 + 800j, 625 + 700j, 725 + 700j, 750 + 800j)
        self.assertAlmostEqual(path8.point(0), (600 + 800j))
        self.assertAlmostEqual(path8.point(0.3), (638.7 + 737j))
        self.assertAlmostEqual(path8.point(0.5), (675 + 725j))
        self.assertAlmostEqual(path8.point(0.9), (740.4 + 773j))
        self.assertAlmostEqual(path8.point(1), (750 + 800j))

        # S875,900 900,800
        inversion = (750 + 800j) + (750 + 800j) - (725 + 700j)
        path9 = CubicBezier(750 + 800j, inversion, 875 + 900j, 900 + 800j)
        self.assertAlmostEqual(path9.point(0), (750 + 800j))
        self.assertAlmostEqual(path9.point(0.3), (788.7 + 863j))
        self.assertAlmostEqual(path9.point(0.5), (825 + 875j))
        self.assertAlmostEqual(path9.point(0.9), (890.4 + 827j))
        self.assertAlmostEqual(path9.point(1), (900 + 800j))

    def test_length(self):

        # A straight line:
        cub = CubicBezier(
            complex(0, 0),
            complex(0, 0),
            complex(0, 100),
            complex(0, 100)
        )

        self.assertAlmostEqual(cub.length(), 100)

        # A diagonal line:
        cub = CubicBezier(
            complex(0, 0),
            complex(0, 0),
            complex(100, 100),
            complex(100, 100)
        )

        self.assertAlmostEqual(cub.length(), sqrt(2 * 100 * 100))

        # A quarter circle large_arc with radius 100:
        kappa = 4 * (sqrt(2) - 1) / 3  # http://www.whizkidtech.redprince.net/bezier/circle/

        cub = CubicBezier(
            complex(0, 0),
            complex(0, kappa * 100),
            complex(100 - kappa * 100, 100),
            complex(100, 100)
        )

        # We can't compare with pi*50 here, because this is just an
        # approximation of a circle large_arc. pi*50 is 157.079632679
        # So this is just yet another "warn if this changes" test.
        # This value is not verified to be correct.
        self.assertAlmostEqual(cub.length(), 157.1016698)

        # A recursive solution has also been suggested, but for CubicBezier
        # curves it could get a false solution on curves where the midpoint is
        # on a straight line between the start and end. For example, the
        # following curve would get solved as a straight line and get the
        # length 300.
        # Make sure this is not the case.
        cub = CubicBezier(
            complex(600, 500),
            complex(600, 350),
            complex(900, 650),
            complex(900, 500)
        )
        self.assertTrue(cub.length() > 300.0)

    def test_equality(self):
        # This is to test the __eq__ and __ne__ methods, so we can't use
        # assertEqual and assertNotEqual
        segment = CubicBezier(complex(600, 500), complex(600, 350),
                              complex(900, 650), complex(900, 500))

        self.assertTrue(segment ==
                        CubicBezier(600 + 500j, 600 + 350j, 900 + 650j, 900 + 500j))
        self.assertTrue(segment !=
                        CubicBezier(600 + 501j, 600 + 350j, 900 + 650j, 900 + 500j))
        self.assertTrue(segment != Line(0, 400))


class QuadraticBezierTest(unittest.TestCase):

    def test_svg_examples(self):
        """These is the path in the SVG specs"""
        # M200,300 Q400,50 600,300 T1000,300
        path1 = QuadraticBezier(200 + 300j, 400 + 50j, 600 + 300j)
        self.assertAlmostEqual(path1.point(0), (200 + 300j))
        self.assertAlmostEqual(path1.point(0.3), (320 + 195j))
        self.assertAlmostEqual(path1.point(0.5), (400 + 175j))
        self.assertAlmostEqual(path1.point(0.9), (560 + 255j))
        self.assertAlmostEqual(path1.point(1), (600 + 300j))

        # T1000, 300
        inversion = (600 + 300j) + (600 + 300j) - (400 + 50j)
        path2 = QuadraticBezier(600 + 300j, inversion, 1000 + 300j)
        self.assertAlmostEqual(path2.point(0), (600 + 300j))
        self.assertAlmostEqual(path2.point(0.3), (720 + 405j))
        self.assertAlmostEqual(path2.point(0.5), (800 + 425j))
        self.assertAlmostEqual(path2.point(0.9), (960 + 345j))
        self.assertAlmostEqual(path2.point(1), (1000 + 300j))

    def test_length(self):
        # expected results calculated with
        # svg.path.segment_length(q, 0, 1, q.start, q.end, 1e-14, 20, 0)
        q1 = QuadraticBezier(200 + 300j, 400 + 50j, 600 + 300j)
        q2 = QuadraticBezier(200 + 300j, 400 + 50j, 500 + 200j)
        closedq = QuadraticBezier(6+2j, 5-1j, 6+2j)
        linq1 = QuadraticBezier(1, 2, 3)
        linq2 = QuadraticBezier(1+3j, 2+5j, -9 - 17j)
        nodalq = QuadraticBezier(1, 1, 1)
        tests = [(q1, 487.77109389525975),
                 (q2, 379.90458193489155),
                 (closedq, 3.1622776601683795),
                 (linq1, 2),
                 (linq2, 22.73335777124786),
                 (nodalq, 0)]
        for q, exp_res in tests:
            self.assertAlmostEqual(q.length(), exp_res)

        # partial length tests
        tests = [(q1, 212.34775387566032),
                 (q2, 166.22170622052397),
                 (closedq, 0.7905694150420949),
                 (linq1, 1.0),
                 (nodalq, 0.0)]
        t0 = 0.25
        t1 = 0.75
        for q, exp_res in tests:
            self.assertAlmostEqual(q.length(t0=t0, t1=t1), exp_res)

        # linear partial cases
        linq2 = QuadraticBezier(1+3j, 2+5j, -9 - 17j)
        tests = [(0, 1/24, 0.13975424859373725),
                 (0, 1/12, 0.1863389981249823),
                 (0, 0.5, 4.844813951249543),
                 (0, 1, 22.73335777124786),
                 (1/24, 1/12, 0.04658474953124506),
                 (1/24, 0.5, 4.705059702655722),
                 (1/24, 1, 22.59360352265412),
                 (1/12, 0.5, 4.658474953124562),
                 (1/12, 1, 22.54701877312288),
                 (0.5, 1, 17.88854381999832)]
        for t0, t1, exp_s in tests:
            self.assertAlmostEqual(linq2.length(t0=t0, t1=t1), exp_s)

    def test_equality(self):
        # This is to test the __eq__ and __ne__ methods, so we can't use
        # assertEqual and assertNotEqual
        segment = QuadraticBezier(200 + 300j, 400 + 50j, 600 + 300j)
        self.assertTrue(segment == QuadraticBezier(200 + 300j, 400 + 50j, 600 + 300j))
        self.assertTrue(segment != QuadraticBezier(200 + 301j, 400 + 50j, 600 + 300j))
        self.assertFalse(segment == Arc(0j, 100 + 50j, 0, 0, 0, 100 + 50j))
        self.assertTrue(Arc(0j, 100 + 50j, 0, 0, 0, 100 + 50j) != segment)


class ArcTest(unittest.TestCase):

    def test_trusting_acos(self):
        """`u1.real` is > 1 in this arc due to numerical error."""
        try:
            a1 = Arc(start=(160.197+102.925j), 
                     radius=(0.025+0.025j), 
                     rotation=0.0, 
                     large_arc=False, 
                     sweep=True, 
                     end=(160.172+102.95j))
        except ValueError:
            self.fail("Arc() raised ValueError unexpectedly!")

    def test_points(self):
        arc1 = Arc(0j, 100 + 50j, 0, 0, 0, 100 + 50j)
        self.assertAlmostEqual(arc1.center, 100 + 0j)
        self.assertAlmostEqual(arc1.theta, 180.0)
        self.assertAlmostEqual(arc1.delta, -90.0)

        self.assertAlmostEqual(arc1.point(0.0), 0j)
        self.assertAlmostEqual(arc1.point(0.1), (1.23116594049 + 7.82172325201j))
        self.assertAlmostEqual(arc1.point(0.2), (4.89434837048 + 15.4508497187j))
        self.assertAlmostEqual(arc1.point(0.3), (10.8993475812 + 22.699524987j))
        self.assertAlmostEqual(arc1.point(0.4), (19.0983005625 + 29.3892626146j))
        self.assertAlmostEqual(arc1.point(0.5), (29.2893218813 + 35.3553390593j))
        self.assertAlmostEqual(arc1.point(0.6), (41.2214747708 + 40.4508497187j))
        self.assertAlmostEqual(arc1.point(0.7), (54.6009500260 + 44.5503262094j))
        self.assertAlmostEqual(arc1.point(0.8), (69.0983005625 + 47.5528258148j))
        self.assertAlmostEqual(arc1.point(0.9), (84.3565534960 + 49.3844170298j))
        self.assertAlmostEqual(arc1.point(1.0), (100 + 50j))

        arc2 = Arc(0j, 100 + 50j, 0, 1, 0, 100 + 50j)
        self.assertAlmostEqual(arc2.center, 50j)
        self.assertAlmostEqual(arc2.theta, -90.0)
        self.assertAlmostEqual(arc2.delta, -270.0)

        self.assertAlmostEqual(arc2.point(0.0), 0j)
        self.assertAlmostEqual(arc2.point(0.1), (-45.399049974 + 5.44967379058j))
        self.assertAlmostEqual(arc2.point(0.2), (-80.9016994375 + 20.6107373854j))
        self.assertAlmostEqual(arc2.point(0.3), (-98.7688340595 + 42.178276748j))
        self.assertAlmostEqual(arc2.point(0.4), (-95.1056516295 + 65.4508497187j))
        self.assertAlmostEqual(arc2.point(0.5), (-70.7106781187 + 85.3553390593j))
        self.assertAlmostEqual(arc2.point(0.6), (-30.9016994375 + 97.5528258148j))
        self.assertAlmostEqual(arc2.point(0.7), (15.643446504 + 99.3844170298j))
        self.assertAlmostEqual(arc2.point(0.8), (58.7785252292 + 90.4508497187j))
        self.assertAlmostEqual(arc2.point(0.9), (89.1006524188 + 72.699524987j))
        self.assertAlmostEqual(arc2.point(1.0), (100 + 50j))

        arc3 = Arc(0j, 100 + 50j, 0, 0, 1, 100 + 50j)
        self.assertAlmostEqual(arc3.center, 50j)
        self.assertAlmostEqual(arc3.theta, -90.0)
        self.assertAlmostEqual(arc3.delta, 90.0)

        self.assertAlmostEqual(arc3.point(0.0), 0j)
        self.assertAlmostEqual(arc3.point(0.1), (15.643446504 + 0.615582970243j))
        self.assertAlmostEqual(arc3.point(0.2), (30.9016994375 + 2.44717418524j))
        self.assertAlmostEqual(arc3.point(0.3), (45.399049974 + 5.44967379058j))
        self.assertAlmostEqual(arc3.point(0.4), (58.7785252292 + 9.54915028125j))
        self.assertAlmostEqual(arc3.point(0.5), (70.7106781187 + 14.6446609407j))
        self.assertAlmostEqual(arc3.point(0.6), (80.9016994375 + 20.6107373854j))
        self.assertAlmostEqual(arc3.point(0.7), (89.1006524188 + 27.300475013j))
        self.assertAlmostEqual(arc3.point(0.8), (95.1056516295 + 34.5491502813j))
        self.assertAlmostEqual(arc3.point(0.9), (98.7688340595 + 42.178276748j))
        self.assertAlmostEqual(arc3.point(1.0), (100 + 50j))

        arc4 = Arc(0j, 100 + 50j, 0, 1, 1, 100 + 50j)
        self.assertAlmostEqual(arc4.center, 100 + 0j)
        self.assertAlmostEqual(arc4.theta, 180.0)
        self.assertAlmostEqual(arc4.delta, 270.0)

        self.assertAlmostEqual(arc4.point(0.0), 0j)
        self.assertAlmostEqual(arc4.point(0.1), (10.8993475812 - 22.699524987j))
        self.assertAlmostEqual(arc4.point(0.2), (41.2214747708 - 40.4508497187j))
        self.assertAlmostEqual(arc4.point(0.3), (84.3565534960 - 49.3844170298j))
        self.assertAlmostEqual(arc4.point(0.4), (130.901699437 - 47.5528258148j))
        self.assertAlmostEqual(arc4.point(0.5), (170.710678119 - 35.3553390593j))
        self.assertAlmostEqual(arc4.point(0.6), (195.105651630 - 15.4508497187j))
        self.assertAlmostEqual(arc4.point(0.7), (198.768834060 + 7.82172325201j))
        self.assertAlmostEqual(arc4.point(0.8), (180.901699437 + 29.3892626146j))
        self.assertAlmostEqual(arc4.point(0.9), (145.399049974 + 44.5503262094j))
        self.assertAlmostEqual(arc4.point(1.0), (100 + 50j))

        arc5 = Arc((725.307482225571-915.5548199281527j),
                   (202.79421639137703+148.77294617167183j),
                   225.6910319606926, 1, 1,
                   (-624.6375539637027+896.5483089399895j))
        self.assertAlmostEqual(arc5.point(0.0), (725.307482226-915.554819928j))
        self.assertAlmostEqual(arc5.point(0.0909090909091), (1023.47397369-597.730444283j))
        self.assertAlmostEqual(arc5.point(0.181818181818), (1242.80253007-232.251400124j))
        self.assertAlmostEqual(arc5.point(0.272727272727), (1365.52445614+151.273373978j))
        self.assertAlmostEqual(arc5.point(0.363636363636), (1381.69755131+521.772981736j))
        self.assertAlmostEqual(arc5.point(0.454545454545), (1290.01156757+849.231748376j))
        self.assertAlmostEqual(arc5.point(0.545454545455), (1097.89435807+1107.12091209j))
        self.assertAlmostEqual(arc5.point(0.636363636364), (820.910116547+1274.54782658j))
        self.assertAlmostEqual(arc5.point(0.727272727273), (481.49845896+1337.94855893j))
        self.assertAlmostEqual(arc5.point(0.818181818182), (107.156499251+1292.18675889j))
        self.assertAlmostEqual(arc5.point(0.909090909091), (-271.788803303+1140.96977533j))

    def test_length(self):
        # I'll test the length calculations by making a circle, in two parts.
        arc1 = Arc(0j, 100 + 100j, 0, 0, 0, 200 + 0j)
        arc2 = Arc(200 + 0j, 100 + 100j, 0, 0, 0, 0j)
        self.assertAlmostEqual(arc1.length(), pi * 100)
        self.assertAlmostEqual(arc2.length(), pi * 100)

    def test_equality(self):
        # This is to test the __eq__ and __ne__ methods, so we can't use
        # assertEqual and assertNotEqual
        segment = Arc(0j, 100 + 50j, 0, 0, 0, 100 + 50j)
        self.assertTrue(segment == Arc(0j, 100 + 50j, 0, 0, 0, 100 + 50j))
        self.assertTrue(segment != Arc(0j, 100 + 50j, 0, 1, 0, 100 + 50j))

    def test_point_to_t(self):
        a = Arc(start=(0+0j), radius=(5+5j), rotation=0.0, large_arc=True, sweep=True, end=(0+10j))
        assert(a.point_to_t(0+0j) == 0.0)
        assert(isclose(a.point_to_t(5+5j), 0.5))
        assert(a.point_to_t(0+10j) == 1.0)
        assert(a.point_to_t(-5+5j) == None)
        assert(a.point_to_t(0+5j) == None)
        assert(a.point_to_t(1+0j) == None)
        assert(a.point_to_t(0-1j) == None)
        assert(a.point_to_t(0+11j) == None)

        a = Arc(start=(0+0j), radius=(5+5j), rotation=0.0, large_arc=True, sweep=False, end=(0+10j))
        assert(a.point_to_t(0+0j) == 0.0)
        assert(isclose(a.point_to_t(-5+5j), 0.5))
        assert(a.point_to_t(0+10j) == 1.0)
        assert(a.point_to_t(5+5j) == None)
        assert(a.point_to_t(0+5j) == None)
        assert(a.point_to_t(1+0j) == None)
        assert(a.point_to_t(0-1j) == None)
        assert(a.point_to_t(0+11j) == None)

        a = Arc(start=(-10+0j), radius=(10+20j), rotation=0.0, large_arc=True, sweep=True, end=(10+0j))
        assert(a.point_to_t(-10+0j) == 0.0)
        assert(isclose(a.point_to_t(0-20j), 0.5))
        assert(a.point_to_t(10+0j) == 1.0)
        assert(a.point_to_t(0+20j) == None)
        assert(a.point_to_t(-5+5j) == None)
        assert(a.point_to_t(0+5j) == None)
        assert(a.point_to_t(1+0j) == None)
        assert(a.point_to_t(0-1j) == None)
        assert(a.point_to_t(0+11j) == None)

        a = Arc(start=(100.834+27.987j), radius=(60.6+60.6j), rotation=0.0, large_arc=False, sweep=False, end=(40.234-32.613j))
        assert(a.point_to_t(100.834+27.987j) == 0.0)
        assert(isclose(a.point_to_t(96.2210993246+4.7963831644j), 0.25))
        assert(isclose(a.point_to_t(83.0846703014-14.8636715784j), 0.50))
        assert(isclose(a.point_to_t(63.4246151671-28.0001000158j), 0.75))
        assert(a.point_to_t(40.234-32.613j) == 1.00)
        assert(a.point_to_t(-10+0j) == None)
        assert(a.point_to_t(0+0j) == None)

        a = Arc(start=(423.049961698-41.3779390229j), radius=(904.283878032+597.298520765j), rotation=0.0, large_arc=True, sweep=False, end=(548.984030235-312.385118044j))
        orig_t = 0.854049465076
        p = a.point(orig_t)
        computed_t = a.point_to_t(p)
        assert(isclose(orig_t, computed_t))

        a = Arc(start=(-1-750j), radius=(750+750j), rotation=0.0, large_arc=True, sweep=False, end=1-750j)
        assert(isclose(a.point_to_t(730.5212132777968+169.8191111892562j), 0.71373858))
        assert(a.point_to_t(730.5212132777968+169j) == None)
        assert(a.point_to_t(730.5212132777968+171j) == None)

        random.seed()
        for arc_index in range(100):
            a = random_arc()
            print(a)
            for t_index in range(100):
                orig_t = random.random()
                p = a.point(orig_t)
                computed_t = a.point_to_t(p)
                print("t:", orig_t)
                print("p:", p)
                print("computed t:", computed_t)
                assert(isclose(orig_t, computed_t))


class TestPath(unittest.TestCase):

    def test_circle(self):
        arc1 = Arc(0j, 100 + 100j, 0, 0, 0, 200 + 0j)
        arc2 = Arc(200 + 0j, 100 + 100j, 0, 0, 0, 0j)
        path = Path(arc1, arc2)
        self.assertAlmostEqual(path.point(0.0), 0j)
        self.assertAlmostEqual(path.point(0.25), (100 + 100j))
        self.assertAlmostEqual(path.point(0.5), (200 + 0j))
        self.assertAlmostEqual(path.point(0.75), (100 - 100j))
        self.assertAlmostEqual(path.point(1.0), 0j)
        self.assertAlmostEqual(path.length(), pi * 200)

    def test_svg_specs(self):
        """The paths that are in the SVG specs"""

        # Big pie: M300,200 h-150 a150,150 0 1,0 150,-150 z
        path = Path(Line(300 + 200j, 150 + 200j),
                    Arc(150 + 200j, 150 + 150j, 0, 1, 0, 300 + 50j),
                    Line(300 + 50j, 300 + 200j))
        # The points and length for this path are calculated and not regression tests.
        self.assertAlmostEqual(path.point(0.0), (300 + 200j))
        self.assertAlmostEqual(path.point(0.14897825542), (150 + 200j))
        self.assertAlmostEqual(path.point(0.5), (406.066017177 + 306.066017177j))
        self.assertAlmostEqual(path.point(1 - 0.14897825542), (300 + 50j))
        self.assertAlmostEqual(path.point(1.0), (300 + 200j))
        # The errors seem to accumulate. Still 6 decimal places is more than good enough.
        self.assertAlmostEqual(path.length(), pi * 225 + 300, places=6)

        # Little pie: M275,175 v-150 a150,150 0 0,0 -150,150 z
        path = Path(Line(275 + 175j, 275 + 25j),
                    Arc(275 + 25j, 150 + 150j, 0, 0, 0, 125 + 175j),
                    Line(125 + 175j, 275 + 175j))
        # The points and length for this path are calculated and not regression tests.
        self.assertAlmostEqual(path.point(0.0), (275 + 175j))
        self.assertAlmostEqual(path.point(0.2800495767557787), (275 + 25j))
        self.assertAlmostEqual(path.point(0.5), (168.93398282201787 + 68.93398282201787j))
        self.assertAlmostEqual(path.point(1 - 0.2800495767557787), (125 + 175j))
        self.assertAlmostEqual(path.point(1.0), (275 + 175j))
        # The errors seem to accumulate. Still 6 decimal places is more than good enough.
        self.assertAlmostEqual(path.length(), pi * 75 + 300, places=6)

        # Bumpy path: M600,350 l 50,-25
        #             a25,25 -30 0,1 50,-25 l 50,-25
        #             a25,50 -30 0,1 50,-25 l 50,-25
        #             a25,75 -30 0,1 50,-25 l 50,-25
        #             a25,100 -30 0,1 50,-25 l 50,-25

        # Commented out because by Andy cause I was skeptical of path.point
        # ground truth values
        # path = Path(Line(600 + 350j, 650 + 325j),
        #             Arc(650 + 325j, 25 + 25j, -30, 0, 1, 700 + 300j),
        #             Line(700 + 300j, 750 + 275j),
        #             Arc(750 + 275j, 25 + 50j, -30, 0, 1, 800 + 250j),
        #             Line(800 + 250j, 850 + 225j),
        #             Arc(850 + 225j, 25 + 75j, -30, 0, 1, 900 + 200j),
        #             Line(900 + 200j, 950 + 175j),
        #             Arc(950 + 175j, 25 + 100j, -30, 0, 1, 1000 + 150j),
        #             Line(1000 + 150j, 1050 + 125j),
        #             )
        # # These are *not* calculated, but just regression tests. Be skeptical.
        # self.assertAlmostEqual(path.point(0), (600+350j))
        # self.assertAlmostEqual(path.point(0.3), (755.239799276+212.182020958j))
        # self.assertAlmostEqual(path.point(0.5), (827.730749264+147.824157418j))
        # self.assertAlmostEqual(path.point(0.9), (971.284357806+106.302352605j))
        # self.assertAlmostEqual(path.point(1), (1050+125j))
        # # The errors seem to accumulate. Still 6 decimal places is more than good enough.
        # self.assertAlmostEqual(path.length(), 928.3886394081095)

    def test_repr(self):
        path = Path(
            Line(start=600 + 350j, end=650 + 325j),
            Arc(start=650 + 325j, radius=25 + 25j, rotation=-30, large_arc=0, sweep=1, end=700 + 300j),
            CubicBezier(start=700 + 300j, control1=800 + 400j, control2=750 + 200j, end=600 + 100j),
            QuadraticBezier(start=600 + 100j, control=600, end=600 + 300j))
        self.assertEqual(eval(repr(path)), path)

    def test_equality(self):
        # This is to test the __eq__ and __ne__ methods, so we can't use
        # assertEqual and assertNotEqual
        path1 = Path(
            Line(start=600 + 350j, end=650 + 325j),
            Arc(start=650 + 325j, radius=25 + 25j, rotation=-30, large_arc=0, sweep=1, end=700 + 300j),
            CubicBezier(start=700 + 300j, control1=800 + 400j, control2=750 + 200j, end=600 + 100j),
            QuadraticBezier(start=600 + 100j, control=600, end=600 + 300j))
        path2 = Path(
            Line(start=600 + 350j, end=650 + 325j),
            Arc(start=650 + 325j, radius=25 + 25j, rotation=-30, large_arc=0, sweep=1, end=700 + 300j),
            CubicBezier(start=700 + 300j, control1=800 + 400j, control2=750 + 200j, end=600 + 100j),
            QuadraticBezier(start=600 + 100j, control=600, end=600 + 300j))

        self.assertTrue(path1 == path2)
        # Modify path2:
        path2[0].start = 601 + 350j
        self.assertTrue(path1 != path2)

        # Modify back:
        path2[0].start = 600 + 350j
        self.assertFalse(path1 != path2)

        # Get rid of the last segment:
        del path2[-1]
        self.assertFalse(path1 == path2)

        # It's not equal to a list of it's segments
        self.assertTrue(path1 != path1[:])
        self.assertFalse(path1 == path1[:])

    def test_continuous_subpaths(self):
        """Test the Path.continuous_subpaths() method."""

        # Continuous and open example
        q = Path(Line(1, 2))
        a = [Path(Line(1, 2))]
        subpaths = q.continuous_subpaths()
        chk1 = all(subpath.iscontinuous() for subpath in subpaths)
        chk2 = (q == Path(*[seg for subpath in subpaths for seg in subpath]))
        self.assertTrue(subpaths == a)
        self.assertTrue(chk1)
        self.assertTrue(chk2)

        # # Continuous and closed example
        q = Path(Line(1, 2), Line(2, 1))
        a = [Path(Line(1, 2), Line(2, 1))]
        subpaths = q.continuous_subpaths()
        chk1 = all(subpath.iscontinuous() for subpath in subpaths)
        chk2 = q == Path(*[seg for subpath in subpaths for seg in subpath])
        self.assertTrue(subpaths == a)
        self.assertTrue(chk1)
        self.assertTrue(chk2)

        # Continuous and open example
        q = Path(Line(1, 2), Line(2, 3), Line(3, 4))
        a = [Path(Line(1, 2), Line(2, 3), Line(3, 4))]
        subpaths = q.continuous_subpaths()
        chk1 = all(subpath.iscontinuous() for subpath in subpaths)
        chk2 = (q == Path(*[seg for subpath in subpaths for seg in subpath]))
        self.assertTrue(subpaths == a)
        self.assertTrue(chk1)
        self.assertTrue(chk2)

        # Continuous and closed example
        q = Path(Line(1, 2), Line(2, 3), Line(3, 4), Line(4, 1))
        a = [Path(Line(1, 2), Line(2, 3), Line(3, 4), Line(4, 1))]
        subpaths = q.continuous_subpaths()
        chk1 = all(subpath.iscontinuous() for subpath in subpaths)
        chk2 = (q == Path(*[seg for subpath in subpaths for seg in subpath]))
        self.assertTrue(subpaths == a)
        self.assertTrue(chk1)
        self.assertTrue(chk2)

        # Discontinuous example
        q = Path(Line(1, 2), Line(2, 3), Line(3, 4),
                 Line(10, 11))
        a = [Path(Line(1, 2), Line(2, 3), Line(3, 4)),
             Path(Line(10, 11))]
        subpaths = q.continuous_subpaths()
        chk1 = all(subpath.iscontinuous() for subpath in subpaths)
        chk2 = (q == Path(*[seg for subpath in subpaths for seg in subpath]))
        self.assertTrue(subpaths == a)
        self.assertTrue(chk1)
        self.assertTrue(chk2)

        # Discontinuous closed example
        q = Path(Line(1, 2), Line(2, 3), Line(3, 4), Line(4, 1),
                 Line(10, 11), Line(11, 12))
        a = [Path(Line(1, 2), Line(2, 3), Line(3, 4), Line(4, 1)),
             Path(Line(10, 11), Line(11, 12))]
        subpaths = q.continuous_subpaths()
        chk1 = all(subpath.iscontinuous() for subpath in subpaths)
        chk2 = (q == Path(*[seg for subpath in subpaths for seg in subpath]))
        self.assertTrue(subpaths == a)
        self.assertTrue(chk1)
        self.assertTrue(chk2)

        # Discontinuous example
        q = Path(Line(1, 2),
                 Line(1, 2), Line(2, 3),
                 Line(10, 11), Line(11, 12), Line(12, 13),
                 Line(10, 11), Line(11, 12), Line(12, 13), Line(13, 14))
        a = [Path(Line(1, 2)),
             Path(Line(1, 2), Line(2, 3)),
             Path(Line(10, 11), Line(11, 12), Line(12, 13)),
             Path(Line(10, 11), Line(11, 12), Line(12, 13), Line(13, 14))]
        subpaths = q.continuous_subpaths()
        chk1 = all(subpath.iscontinuous() for subpath in subpaths)
        chk2 = (q == Path(*[seg for subpath in subpaths for seg in subpath]))
        self.assertTrue(subpaths == a)
        self.assertTrue(chk1)
        self.assertTrue(chk2)

        # Discontinuous example with overlapping end
        q = Path(Line(1, 2),
                 Line(5, 6), Line(6, 7),
                 Line(10, 11), Line(11, 12), Line(12, 13),
                 Line(10, 11), Line(11, 12), Line(12, 13), Line(13, 1))
        a = [Path(Line(1, 2)),
             Path(Line(5, 6), Line(6, 7)),
             Path(Line(10, 11), Line(11, 12), Line(12, 13)),
             Path(Line(10, 11), Line(11, 12), Line(12, 13), Line(13, 1))]
        subpaths = q.continuous_subpaths()
        chk1 = all(subpath.iscontinuous() for subpath in subpaths)
        chk2 = (q == Path(*[seg for subpath in subpaths for seg in subpath]))
        self.assertTrue(subpaths == a)
        self.assertTrue(chk1)
        self.assertTrue(chk2)

    def test_cropped(self):
        p_closed = Path(Line(0, 1), Line(1, 1 + 1j), Line(1 + 1j, 1j),
                        Line(1j, 0))
        first_half = Path(Line(0, 1), Line(1, 1 + 1j))
        second_half = Path(Line(1 + 1j, 1j), Line(1j, 0))
        middle_half = Path(Line(1, 1 + 1j), Line(1 + 1j, 1j))
        other_middle_half = Path(Line(1j, 0), Line(0, 1))
        self.assertTrue(p_closed.cropped(0, 0.5) == first_half)
        self.assertTrue(p_closed.cropped(1, 0.5) == first_half)
        self.assertTrue(p_closed.cropped(.5, 1) == second_half)
        self.assertTrue(p_closed.cropped(0.25, 0.75) == middle_half)
        self.assertTrue(p_closed.cropped(0.75, 0.25) == other_middle_half)
        with self.assertRaises(AssertionError):
            p_closed.cropped(1, 0)
        with self.assertRaises(AssertionError):
            p_closed.cropped(.5, 1.1)
        with self.assertRaises(AssertionError):
            p_closed.cropped(-0.1, 0.1)

        p_open = Path(Line(0, 1), Line(1, 1 + 1j), Line(1 + 1j, 1j),
                      Line(1j, 2j))
        self.assertTrue(p_open.cropped(0, 0.5) == first_half)
        with self.assertRaises(ValueError):
            p_open.cropped(.75, .25)
        with self.assertRaises(ValueError):
            p_open.cropped(1, .25)
        with self.assertRaises(AssertionError):
            p_open.cropped(1, 0)


class Test_ilength(unittest.TestCase):
    # See svgpathtools.notes.inv_arclength.py for information on how these
    # test values were generated (using the .length() method).
    ##############################################################

    def test_ilength_lines(self):
        l = Line(1, 3-1j)
        nodall = Line(1+1j, 1+1j)

        tests = [(l, 0.01, 0.022360679774997897),
         (l, 0.1, 0.223606797749979),
         (l, 0.5, 1.118033988749895),
         (l, 0.9, 2.012461179749811),
         (l, 0.99, 2.213707297724792)]

        for (l, t, s) in tests:
            self.assertAlmostEqual(l.ilength(s), t)

    def test_ilength_quadratics(self):
        q1 = QuadraticBezier(200 + 300j, 400 + 50j, 600 + 300j)
        q2 = QuadraticBezier(200 + 300j, 400 + 50j, 500 + 200j)
        closedq = QuadraticBezier(6 + 2j, 5 - 1j, 6 + 2j)
        linq = QuadraticBezier(1+3j, 2+5j, -9 - 17j)
        nodalq = QuadraticBezier(1, 1, 1)

        tests = [(q1, 0.01, 6.364183310105577),
         (q1, 0.1, 60.23857499635088),
         (q1, 0.5, 243.8855469477619),
         (q1, 0.9, 427.53251889917294),
         (q1, 0.99, 481.40691058541813),
         (q2, 0.01, 6.365673533661836),
         (q2, 0.1, 60.31675895732397),
         (q2, 0.5, 233.24592830045907),
         (q2, 0.9, 346.42891253298706),
         (q2, 0.99, 376.32659156736844),
         (closedq, 0.01, 0.06261309767133393),
         (closedq, 0.1, 0.5692099788303084),
         (closedq, 0.5, 1.5811388300841898),
         (closedq, 0.9, 2.5930676813380713),
         (closedq, 0.99, 3.0996645624970456),
         (linq, 0.01, 0.04203807797699605),
         (linq, 0.1, 0.19379255804998186),
         (linq, 0.5, 4.844813951249544),
         (linq, 0.9, 18.0823363780483),
         (linq, 0.99, 22.24410609777091)]

        for q, t, s in tests:
            try:
                self.assertAlmostEqual(q.ilength(s), t)
            except:
                print(q)
                print(s)
                print(t)
                raise

    def test_ilength_cubics(self):
        c1 = CubicBezier(200 + 300j, 400 + 50j, 600+100j, -200)
        symc = CubicBezier(1-2j, 10-1j, 10+1j, 1+2j)
        closedc = CubicBezier(1-2j, 10-1j, 10+1j, 1-2j)

        tests = [(c1, 0.01, 9.53434737943073),
                 (c1, 0.1, 88.89941848775852),
                 (c1, 0.5, 278.5750942713189),
                 (c1, 0.9, 651.4957786584646),
                 (c1, 0.99, 840.2010603832538),
                 (symc, 0.01, 0.2690118556702902),
                 (symc, 0.1, 2.45230693868727),
                 (symc, 0.5, 7.256147083644424),
                 (symc, 0.9, 12.059987228602886),
                 (symc, 0.99, 14.243282311619401),
                 (closedc, 0.01, 0.26901140075538765),
                 (closedc, 0.1, 2.451722765460998),
                 (closedc, 0.5, 6.974058969750422),
                 (closedc, 0.9, 11.41781741489913),
                 (closedc, 0.99, 13.681324783697782)]

        for (c, t, s) in tests:
            self.assertAlmostEqual(c.ilength(s), t)

    def test_ilength_arcs(self):
        arc1 = Arc(0j, 100 + 50j, 0, 0, 0, 100 + 50j)
        arc2 = Arc(0j, 100 + 50j, 0, 1, 0, 100 + 50j)
        arc3 = Arc(0j, 100 + 50j, 0, 0, 1, 100 + 50j)
        arc4 = Arc(0j, 100 + 50j, 0, 1, 1, 100 + 50j)
        arc5 = Arc(0j, 100 + 100j, 0, 0, 0, 200 + 0j)
        arc6 = Arc(200 + 0j, 100 + 100j, 0, 0, 0, 0j)
        arc7 = Arc(0j, 100 + 50j, 0, 0, 0, 100 + 50j)

        tests = [(arc1, 0.01, 0.785495042476231),
                 (arc1, 0.1, 7.949362877455911),
                 (arc1, 0.5, 48.28318721111137),
                 (arc1, 0.9, 105.44598206942156),
                 (arc1, 0.99, 119.53485487631241),
                 (arc2, 0.01, 4.71108115728524),
                 (arc2, 0.1, 45.84152747676626),
                 (arc2, 0.5, 169.38878996795734),
                 (arc2, 0.9, 337.44707303579696),
                 (arc2, 0.99, 360.95800139278765),
                 (arc3, 0.01, 1.5707478805335624),
                 (arc3, 0.1, 15.659620687424416),
                 (arc3, 0.5, 72.82241554573457),
                 (arc3, 0.9, 113.15623987939003),
                 (arc3, 0.99, 120.3201077143697),
                 (arc4, 0.01, 2.3588068777503897),
                 (arc4, 0.1, 25.869735234740887),
                 (arc4, 0.5, 193.9280183025816),
                 (arc4, 0.9, 317.4752807937718),
                 (arc4, 0.99, 358.6057271132536),
                 (arc5, 0.01, 3.141592653589793),
                 (arc5, 0.1, 31.415926535897935),
                 (arc5, 0.5, 157.07963267948966),
                 (arc5, 0.9, 282.7433388230814),
                 (arc5, 0.99, 311.01767270538954),
                 (arc6, 0.01, 3.141592653589793),
                 (arc6, 0.1, 31.415926535897928),
                 (arc6, 0.5, 157.07963267948966),
                 (arc6, 0.9, 282.7433388230814),
                 (arc6, 0.99, 311.01767270538954),
                 (arc7, 0.01, 0.785495042476231),
                 (arc7, 0.1, 7.949362877455911),
                 (arc7, 0.5, 48.28318721111137),
                 (arc7, 0.9, 105.44598206942156),
                 (arc7, 0.99, 119.53485487631241)]

        for (c, t, s) in tests:
            self.assertAlmostEqual(c.ilength(s), t)

    def test_ilength_paths(self):
        line1 = Line(600 + 350j, 650 + 325j)
        arc1 = Arc(650 + 325j, 25 + 25j, -30, 0, 1, 700 + 300j)
        cub1 = CubicBezier(650 + 325j, 25 + 25j, -30, 700 + 300j)
        cub2 = CubicBezier(700 + 300j, 800 + 400j, 750 + 200j, 600 + 100j)
        quad3 = QuadraticBezier(600 + 100j, 600, 600 + 300j)
        linez = Line(600 + 300j, 600 + 350j)

        bezpath = Path(line1, cub1, cub2, quad3)
        bezpathz = Path(line1, cub1, cub2, quad3, linez)
        path = Path(line1, arc1, cub2, quad3)
        pathz = Path(line1, arc1, cub2, quad3, linez)
        lpath = Path(linez)
        qpath = Path(quad3)
        cpath = Path(cub1)
        apath = Path(arc1)

        tests = [(bezpath, 0.0, 0.0),
                 (bezpath, 0.1111111111111111, 286.2533595149515),
                 (bezpath, 0.2222222222222222, 503.8620222915423),
                 (bezpath, 0.3333333333333333, 592.6337135346268),
                 (bezpath, 0.4444444444444444, 644.3880677233315),
                 (bezpath, 0.5555555555555556, 835.0384185011363),
                 (bezpath, 0.6666666666666666, 1172.8729938994575),
                 (bezpath, 0.7777777777777778, 1308.6205983178952),
                 (bezpath, 0.8888888888888888, 1532.8473168900994),
                 (bezpath, 1.0, 1758.2427369258733),
                 (bezpathz, 0.0, 0.0),
                 (bezpathz, 0.1111111111111111, 294.15942308605435),
                 (bezpathz, 0.2222222222222222, 512.4295461513882),
                 (bezpathz, 0.3333333333333333, 594.0779370040138),
                 (bezpathz, 0.4444444444444444, 658.7361976564598),
                 (bezpathz, 0.5555555555555556, 874.1674336581542),
                 (bezpathz, 0.6666666666666666, 1204.2371344392693),
                 (bezpathz, 0.7777777777777778, 1356.773042865213),
                 (bezpathz, 0.8888888888888888, 1541.808492602876),
                 (bezpathz, 1.0, 1808.2427369258733),
                 (path, 0.0, 0.0),
                 (path, 0.1111111111111111, 81.44016397108298),
                 (path, 0.2222222222222222, 164.72556816469307),
                 (path, 0.3333333333333333, 206.71343564679154),
                 (path, 0.4444444444444444, 265.4898349999353),
                 (path, 0.5555555555555556, 367.5420981413199),
                 (path, 0.6666666666666666, 487.29863861165995),
                 (path, 0.7777777777777778, 511.84069655405284),
                 (path, 0.8888888888888888, 579.9530841780238),
                 (path, 1.0, 732.9614757397469),
                 (pathz, 0.0, 0.0),
                 (pathz, 0.1111111111111111, 86.99571952663854),
                 (pathz, 0.2222222222222222, 174.33662608180325),
                 (pathz, 0.3333333333333333, 214.42194393858466),
                 (pathz, 0.4444444444444444, 289.94661033436205),
                 (pathz, 0.5555555555555556, 408.38391100702125),
                 (pathz, 0.6666666666666666, 504.4309373835351),
                 (pathz, 0.7777777777777778, 533.774834546298),
                 (pathz, 0.8888888888888888, 652.931321760894),
                 (pathz, 1.0, 782.9614757397469),
                 (lpath, 0.0, 0.0),
                 (lpath, 0.1111111111111111, 5.555555555555555),
                 (lpath, 0.2222222222222222, 11.11111111111111),
                 (lpath, 0.3333333333333333, 16.666666666666664),
                 (lpath, 0.4444444444444444, 22.22222222222222),
                 (lpath, 0.5555555555555556, 27.77777777777778),
                 (lpath, 0.6666666666666666, 33.33333333333333),
                 (lpath, 0.7777777777777778, 38.88888888888889),
                 (lpath, 0.8888888888888888, 44.44444444444444),
                 (lpath, 1.0, 50.0),
                 (qpath, 0.0, 0.0),
                 (qpath, 0.1111111111111111, 17.28395061728395),
                 (qpath, 0.2222222222222222, 24.69135802469136),
                 (qpath, 0.3333333333333333, 27.777777777777786),
                 (qpath, 0.4444444444444444, 40.12345679012344),
                 (qpath, 0.5555555555555556, 62.3456790123457),
                 (qpath, 0.6666666666666666, 94.44444444444446),
                 (qpath, 0.7777777777777778, 136.41975308641975),
                 (qpath, 0.8888888888888888, 188.27160493827154),
                 (qpath, 1.0, 250.0),
                 (cpath, 0.0, 0.0),
                 (cpath, 0.1111111111111111, 207.35525375551356),
                 (cpath, 0.2222222222222222, 366.0583590267552),
                 (cpath, 0.3333333333333333, 474.34064293812787),
                 (cpath, 0.4444444444444444, 530.467036317684),
                 (cpath, 0.5555555555555556, 545.0444351253911),
                 (cpath, 0.6666666666666666, 598.9767847757622),
                 (cpath, 0.7777777777777778, 710.4080903390646),
                 (cpath, 0.8888888888888888, 881.1796899225557),
                 (cpath, 1.0, 1113.0914444911352),
                 (apath, 0.0, 0.0),
                 (apath, 0.1111111111111111, 9.756687033889872),
                 (apath, 0.2222222222222222, 19.51337406777974),
                 (apath, 0.3333333333333333, 29.27006110166961),
                 (apath, 0.4444444444444444, 39.02674813555948),
                 (apath, 0.5555555555555556, 48.783435169449355),
                 (apath, 0.6666666666666666, 58.54012220333922),
                 (apath, 0.7777777777777778, 68.2968092372291),
                 (apath, 0.8888888888888888, 78.05349627111896),
                 (apath, 1.0, 87.81018330500885)]

        for (c, t, s) in tests:
            self.assertAlmostEqual(c.ilength(s), t, msg=str((c, t, s)))

    # Exceptional Cases
    def test_ilength_exceptions(self):
        nodalq = QuadraticBezier(1, 1, 1)
        with self.assertRaises(AssertionError):
            nodalq.ilength(1)

        lin = Line(0, 0.5j)
        with self.assertRaises(ValueError):
            lin.ilength(1)


class Test_intersect(unittest.TestCase):
    def test_intersect(self):

        ###################################################################
        # test that `some_seg.intersect(another_seg)` will produce properly
        # ordered tuples, i.e. the first element in each tuple refers to
        # `some_seg` and the second element refers to `another_seg`.
        # Also tests that the correct number of intersections is found.
        a = Line(0 + 200j, 300 + 200j)
        b = QuadraticBezier(40 + 150j, 70 + 200j, 210 + 300j)
        c = CubicBezier(60 + 150j, 40 + 200j, 120 + 250j, 200 + 160j)
        d = Arc(70 + 150j, 50 + 100j, 0, 0, 0, 200 + 100j)
        segdict = {'line': a, "quadratic": b, 'cubic': c, 'arc': d}

        # test each segment type against each other type
        for x, y in [(x, y) for x in segdict for y in segdict]:
            if x == y:
                continue
            x = segdict[x]
            y = segdict[y]
            xiy = sorted(x.intersect(y, tol=1e-15))
            yix = sorted(y.intersect(x, tol=1e-15), key=itemgetter(1))
            for xy, yx in zip(xiy, yix):
                self.assertAlmostEqual(xy[0], yx[1])
                self.assertAlmostEqual(xy[1], yx[0])
                self.assertAlmostEqual(x.point(xy[0]), y.point(yx[0]))
            self.assertTrue(len(xiy), len(yix))
            self.assertTrue(len(xiy), 2)
            self.assertTrue(len(yix), 2)

        # test each segment against another segment of same type
        for x in segdict:
            if x == 'arc':
                # this is an example of the Arc.intersect method not working
                # in call cases.  See docstring for a note on its
                # incomplete implementation.
                continue
            x = segdict[x]
            y = x.rotated(90).translated(5)
            xiy = sorted(x.intersect(y, tol=1e-15))
            yix = sorted(y.intersect(x, tol=1e-15), key=itemgetter(1))
            for xy, yx in zip(xiy, yix):
                self.assertAlmostEqual(xy[0], yx[1])
                self.assertAlmostEqual(xy[1], yx[0])
                self.assertAlmostEqual(x.point(xy[0]), y.point(yx[0]))
            self.assertTrue(len(xiy), len(yix))
            self.assertTrue(len(xiy), 1)
            self.assertTrue(len(yix), 1)
        ###################################################################

    def test_line_line_0(self):
        l0 = Line(start=(25.389999999999997+99.989999999999995j), end=(25.389999999999997+90.484999999999999j))
        l1 = Line(start=(25.390000000000001+84.114999999999995j), end=(25.389999999999997+74.604202137430320j))
        i = l0.intersect(l1)
        assert(len(i)) == 0

    def test_line_line_1(self):
        l0 = Line(start=(-124.705378549+327.696674827j), end=(12.4926214511+121.261674827j))
        l1 = Line(start=(-12.4926214511+121.261674827j), end=(124.705378549+327.696674827j))
        i = l0.intersect(l1)
        assert(len(i)) == 1
        assert(abs(l0.point(i[0][0])-l1.point(i[0][1])) < 1e-9)


    def test_arc_line(self):
        l = Line(start=(-20+1j), end=(20+1j))
        a = Arc(start=(-10+0), radius=(10+10j), rotation=0.0, large_arc=True, sweep=False, end=(10+0j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 2)

        l = Line(start=(-20-1j), end=(20-1j))
        a = Arc(start=(-10+0), radius=(10+10j), rotation=0.0, large_arc=True, sweep=False, end=(10+0j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 0)

        l = Line(start=(-20+1j), end=(20+1j))
        a = Arc(start=(-10+0), radius=(10+10j), rotation=0.0, large_arc=True, sweep=True, end=(10+0j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 0)

        l = Line(start=(-20-1j), end=(20-1j))
        a = Arc(start=(-10+0), radius=(10+10j), rotation=0.0, large_arc=True, sweep=True, end=(10+0j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 2)

        l = Line(start=(-20+0j), end=(20+0j))
        a = Arc(start=(-10+0), radius=(10+10j), rotation=0.0, large_arc=True, sweep=True, end=(10+0j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 2)

        l = Line(start=(-20+0j), end=(20+0j))
        a = Arc(start=(-10+0), radius=(10+10j), rotation=0.0, large_arc=True, sweep=False, end=(10+0j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 2)

        l = Line(start=(-20+10j), end=(20+10j))
        a = Arc(start=(-10+0), radius=(10+10j), rotation=0.0, large_arc=True, sweep=False, end=(10+0j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 1)

        l = Line(start=(229.226097475-282.403591377j), end=(751.681212592+188.907748894j))
        a = Arc(start=(-1-750j), radius=(750+750j), rotation=0.0, large_arc=True, sweep=False, end=(1-750j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 1)

        # end of arc touches start of horizontal line
        l = Line(start=(40.234-32.613j), end=(12.7-32.613j))
        a = Arc(start=(100.834+27.987j), radius=(60.6+60.6j), rotation=0.0, large_arc=False, sweep=False, end=(40.234-32.613j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 1)

        # vertical line, intersects half-arc once
        l = Line(start=(1-100j), end=(1+100j))
        a = Arc(start=(10.0+0j), radius=(10+10j), rotation=0, large_arc=False, sweep=True, end=(-10.0+0j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 1)

        # vertical line, intersects nearly-full arc twice
        l = Line(start=(1-100j), end=(1+100j))
        a = Arc(start=(0.1-10j), radius=(10+10j), rotation=0, large_arc=True, sweep=True, end=(-0.1-10j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 2)

        # vertical line, start of line touches end of arc
        l = Line(start=(15.4+100j), end=(15.4+90.475j))
        a = Arc(start=(25.4+90j), radius=(10+10j), rotation=0, large_arc=False, sweep=True, end=(15.4+100j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 1)

        l = Line(start=(100-60.913j), end=(40+59j))
        a = Arc(start=(100.834+27.987j), radius=(60.6+60.6j), rotation=0.0, large_arc=False, sweep=False, end=(40.234-32.613j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 1)

        l = Line(start=(128.57143 + 380.93364j), end=(300.00001 + 389.505069j))
        a = Arc(start=(214.28572 + 598.07649j), radius=(85.714287 + 108.57143j), rotation=0.0, large_arc=False, sweep=True, end=(128.57143 + 489.50507j))
        intersections = a.intersect(l)
        assert_intersections(a, l, intersections, 0)

        random.seed()
        for arc_index in range(50):
            a = random_arc()
            print(a)
            for line_index in range(100):
                l = random_line()
                print(l)
                intersections = a.intersect(l)
                assert_intersections(a, l, intersections, None)


    def test_intersect_arc_line_1(self):

        """Verify the return value of intersects() when an Arc ends at
        the starting point of a Line."""

        a = Arc(start=(0+0j), radius=(10+10j), rotation=0, large_arc=False,
                sweep=False, end=(10+10j), autoscale_radius=False)
        l = Line(start=(10+10j), end=(20+10j))

        i = a.intersect(l)
        assert(len(i) == 1)
        assert(i[0][0] == 1.0)
        assert(i[0][1] == 0.0)


    def test_intersect_arc_line_2(self):

        """Verify the return value of intersects() when an Arc is pierced
        once by a Line."""

        a = Arc(start=(0+0j), radius=(10+10j), rotation=0, large_arc=False,
                sweep=False, end=(10+10j), autoscale_radius=False)
        l = Line(start=(0+9j), end=(20+9j))

        i = a.intersect(l)
        assert(len(i) == 1)
        assert(i[0][0] >= 0.0)
        assert(i[0][0] <= 1.0)
        assert(i[0][1] >= 0.0)
        assert(i[0][1] <= 1.0)


    def test_intersect_arc_line_3(self):

        """Verify the return value of intersects() when an Arc misses
        a Line, but the circle that the Arc is part of hits the Line."""

        a = Arc(start=(0+0j), radius=(10+10j), rotation=0, large_arc=False,
                sweep=False, end=(10+10j), autoscale_radius=False)
        l = Line(start=(11+100j), end=(11-100j))

        i = a.intersect(l)
        assert(len(i) == 0)


    def test_intersect_arc_line_disjoint_bboxes(self):
        # The arc is very short, which contributes to the problem here.
        l = Line(start=(125.314540561+144.192926144j), end=(125.798713132+144.510685287j))
        a = Arc(start=(128.26640649+146.908463323j), radius=(2+2j),
                rotation=0, large_arc=False, sweep=True,
                end=(128.26640606+146.90846449j))
        i = l.intersect(a)
        assert(i == [])


class TestPathTools(unittest.TestCase):
    # moved from test_pathtools.py

    def setUp(self):
        self.arc1 = Arc(650+325j, 25+25j, -30.0, False, True, 700+300j)
        self.line1 = Line(0, 100+100j)
        self.quadratic1 = QuadraticBezier(100+100j, 150+150j, 300+200j)
        self.cubic1 = CubicBezier(300+200j, 350+400j, 400+425j, 650+325j)
        self.path_of_all_seg_types = Path(self.line1, self.quadratic1,
                                          self.cubic1, self.arc1)
        self.path_of_bezier_seg_types = Path(self.line1, self.quadratic1,
                                             self.cubic1)

    def test_is_bezier_segment(self):

        # False
        self.assertFalse(is_bezier_segment(self.arc1))
        self.assertFalse(is_bezier_segment(self.path_of_bezier_seg_types))

        # True
        self.assertTrue(is_bezier_segment(self.line1))
        self.assertTrue(is_bezier_segment(self.quadratic1))
        self.assertTrue(is_bezier_segment(self.cubic1))

    def test_is_bezier_path(self):

        # False
        self.assertFalse(is_bezier_path(self.path_of_all_seg_types))
        self.assertFalse(is_bezier_path(self.line1))
        self.assertFalse(is_bezier_path(self.quadratic1))
        self.assertFalse(is_bezier_path(self.cubic1))
        self.assertFalse(is_bezier_path(self.arc1))

        # True
        self.assertTrue(is_bezier_path(self.path_of_bezier_seg_types))
        self.assertTrue(is_bezier_path(Path()))

    def test_polynomial2bezier(self):

        def distfcn(tup1, tup2):
            assert len(tup1) == len(tup2)
            return sum((tup1[i]-tup2[i])**2 for i in range(len(tup1)))**0.5

        # Case: Line
        pcoeffs = [(-1.7-2j), (6+2j)]
        p = poly1d(pcoeffs)
        correct_bpoints = [(6+2j), (4.3+0j)]

        # Input poly1d object
        bez = poly2bez(p)
        bpoints = bez.bpoints()
        self.assertAlmostEqual(distfcn(bpoints, correct_bpoints), 0)

        # Input list of coefficients
        bpoints = poly2bez(pcoeffs, return_bpoints=True)
        self.assertAlmostEqual(distfcn(bpoints, correct_bpoints), 0)

        # Case: Quadratic
        pcoeffs = [(29.5+15.5j), (-31-19j), (7.5+5.5j)]
        p = poly1d(pcoeffs)
        correct_bpoints = [(7.5+5.5j), (-8-4j), (6+2j)]

        # Input poly1d object
        bez = poly2bez(p)
        bpoints = bez.bpoints()
        self.assertAlmostEqual(distfcn(bpoints, correct_bpoints), 0)

        # Input list of coefficients
        bpoints = poly2bez(pcoeffs, return_bpoints=True)
        self.assertAlmostEqual(distfcn(bpoints, correct_bpoints), 0)

        # Case: Cubic
        pcoeffs = [(-18.5-12.5j), (34.5+16.5j), (-18-6j), (6+2j)]
        p = poly1d(pcoeffs)
        correct_bpoints = [(6+2j), 0j, (5.5+3.5j), (4+0j)]

        # Input poly1d object
        bez = poly2bez(p)
        bpoints = bez.bpoints()
        self.assertAlmostEqual(distfcn(bpoints, correct_bpoints), 0)

        # Input list of coefficients object
        bpoints = poly2bez(pcoeffs, return_bpoints=True)
        self.assertAlmostEqual(distfcn(bpoints, correct_bpoints), 0)

    def test_bpoints2bezier(self):
        cubic_bpoints = [(6+2j), 0, (5.5+3.5j), (4+0j)]
        quadratic_bpoints = [(6+2j), 0, (5.5+3.5j)]
        line_bpoints = [(6+2j), 0]
        self.assertTrue(isinstance(bpoints2bezier(cubic_bpoints), CubicBezier))
        self.assertTrue(isinstance(bpoints2bezier(quadratic_bpoints),
                                   QuadraticBezier))
        self.assertTrue(isinstance(bpoints2bezier(line_bpoints), Line))
        self.assertSequenceEqual(bpoints2bezier(cubic_bpoints).bpoints(),
                                 cubic_bpoints)
        self.assertSequenceEqual(bpoints2bezier(quadratic_bpoints).bpoints(),
                                 quadratic_bpoints)
        self.assertSequenceEqual(bpoints2bezier(line_bpoints).bpoints(),
                                 line_bpoints)

    # def test_line2pathd(self):
    #     bpoints = (0+1.5j, 100+10j)
    #     line = Line(*bpoints)
    #
    #     # from Line object
    #     pathd = line2pathd(line)
    #     path = parse_path(pathd)
    #     self.assertTrue(path[0] == line)
    #
    #     # from list of bpoints
    #     pathd = line2pathd(bpoints)
    #     path = parse_path(pathd)
    #     self.assertTrue(path[0] == line)
    #
    # def test_cubic2pathd(self):
    #     bpoints = (0+1.5j, 100+10j, 150-155.3j, 0)
    #     cubic = CubicBezier(*bpoints)
    #
    #     # from Line object
    #     pathd = cubic2pathd(cubic)
    #     path = parse_path(pathd)
    #     self.assertTrue(path[0] == cubic)
    #
    #     # from list of bpoints
    #     pathd = cubic2pathd(bpoints)
    #     path = parse_path(pathd)
    #     self.assertTrue(path[0] == cubic)

    def test_closest_point_in_path(self):
        def distfcn(tup1, tup2):
            assert len(tup1) == len(tup2)
            return sum((tup1[i]-tup2[i])**2 for i in range(len(tup1)))**0.5

        # Note: currently the radiialrange method is not implemented for Arc
        # objects
        # test_path = self.path_of_all_seg_types
        # origin = -123 - 123j
        # expected_result = ???
        # self.assertAlmostEqual(min_radius(origin, test_path),
        # expected_result)

        # generic case (where is_bezier_path(test_path) == True)
        test_path = self.path_of_bezier_seg_types
        pt = 300+300j
        expected_result = (29.382522853493143, 0.17477067969145446, 2)
        result = closest_point_in_path(pt, test_path)
        err = distfcn(expected_result, result)
        self.assertAlmostEqual(err, 0)

        # cubic test with multiple valid solutions
        test_path = Path(CubicBezier(1-2j, 10-1j, 10+1j, 1+2j))
        pt = 3
        expected_results = [(1.7191878932122302, 0.90731678233211366, 0),
                            (1.7191878932122304, 0.092683217667886342, 0)]
        result = closest_point_in_path(pt, test_path)
        err = min(distfcn(e_res, result) for e_res in expected_results)
        self.assertAlmostEqual(err, 0)

    def test_farthest_point_in_path(self):
        def distfcn(tup1, tup2):
            assert len(tup1) == len(tup2)
            return sum((tup1[i]-tup2[i])**2 for i in range(len(tup1)))**0.5

        # Note: currently the radiialrange method is not implemented for Arc
        # objects
        # test_path = self.path_of_all_seg_types
        # origin = -123 - 123j
        # expected_result = ???
        # self.assertAlmostEqual(min_radius(origin, test_path),
        # expected_result)

        # boundary test
        test_path = self.path_of_bezier_seg_types
        pt = 300+300j
        expected_result = (424.26406871192853, 0, 0)
        result = farthest_point_in_path(pt, test_path)
        err = distfcn(expected_result, result)
        self.assertAlmostEqual(err, 0)

        # non-boundary test
        test_path = Path(CubicBezier(1-2j, 10-1j, 10+1j, 1+2j))
        pt = 3
        expected_result = (4.75, 0.5, 0)
        result = farthest_point_in_path(pt, test_path)
        err = distfcn(expected_result, result)
        self.assertAlmostEqual(err, 0)

    def test_path_encloses_pt(self):

        line1 = Line(0, 100+100j)
        quadratic1 = QuadraticBezier(100+100j, 150+150j, 300+200j)
        cubic1 = CubicBezier(300+200j, 350+400j, 400+425j, 650+325j)
        line2 = Line(650+325j, 650+10j)
        line3 = Line(650+10j, 0)
        open_bez_path = Path(line1, quadratic1, cubic1)
        closed_bez_path = Path(line1, quadratic1, cubic1, line2, line3)

        inside_pt = 200+20j
        outside_pt1 = 1000+1000j
        outside_pt2 = 800+800j
        boundary_pt = 50+50j

        # Note: currently the intersect() method is not implemented for Arc
        # objects
        # arc1 = Arc(650+325j, 25+25j, -30.0, False, True, 700+300j)
        # closed_path_with_arc = Path(line1, quadratic1, cubic1, arc1)
        # self.assertTrue(
        #     path_encloses_pt(inside_pt, outside_pt2, closed_path_with_arc))

        # True cases
        self.assertTrue(
            path_encloses_pt(inside_pt, outside_pt2, closed_bez_path))
        self.assertTrue(
            path_encloses_pt(boundary_pt, outside_pt2, closed_bez_path))

        # False cases
        self.assertFalse(
            path_encloses_pt(outside_pt1, outside_pt2, closed_bez_path))

        # Exception Cases
        with self.assertRaises(AssertionError):
            path_encloses_pt(inside_pt, outside_pt2, open_bez_path)

        # Display test paths and points
        # ns2d = [inside_pt, outside_pt1, outside_pt2, boundary_pt]
        # ncolors = ['green', 'red', 'orange', 'purple']
        # disvg(closed_path_with_arc, nodes=ns2d, node_colors=ncolors,
        #       openinbrowser=True)
        # disvg(open_bez_path, nodes=ns2d, node_colors=ncolors,
        #       openinbrowser=True)
        # disvg(closed_bez_path, nodes=ns2d, node_colors=ncolors,
        #       openinbrowser=True)


if __name__ == '__main__':
    unittest.main()
