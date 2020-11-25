#!/usr/bin/env python
# coding:utf-8

from svgpathtools.path import Arc, QuadraticBezier, CubicBezier, Line, Path as sPath
from svgpathtools.parser import parse_path
from svgpathtools.paths2svg import wsvg
from pathlib import Path
import unittest
from numpy import isclose


def visualize_path(path, label, addpoints=None, pathnodes=True, colors=None):
    """Create an SVG file named label.svg containing path."""
    direct = './data/tests_intersects_svg_results'
    directory = Path(direct)
    # it seems that xmldom cannot use PosixPath
    filename = str(directory / str(label + '.svg'))
    if pathnodes:
        startends = [segment.end for segment in path[0:]
                     if not isinstance(segment, Arc)]
        arcs = [segment.end for segment in path[0:]
                if isinstance(segment, Arc)]
        startends.append(path[0].start)
        arcs.append(path[0].start)
        control1s = [segment.control1 for segment in path[0:]
                     if isinstance(segment, CubicBezier)]
        control2s = [segment.control2 for segment in path[0:]
                     if isinstance(segment, CubicBezier)]
        controls = [segment.control for segment in path[0:]
                    if isinstance(segment, QuadraticBezier)]
        controls = control1s + control2s + controls
        nodes = startends + arcs + controls
        node_colors = ['black'] * len(nodes)
    else:
        nodes = []
        node_colors = []
    if addpoints and isinstance(addpoints, complex):
        for point in addpoints:
            nodes.append(point)
        lenaddpoints = len(addpoints)
        node_colors += ['red'] * lenaddpoints
    elif addpoints and isinstance(addpoints, list):
        pointscolors = ['green', 'blue', 'black']
        for i, lis in enumerate(addpoints):
            nodes += lis
            node_colors += [pointscolors[i]] * len(lis)

    wsvg(path,
         filename=filename,
         nodes=nodes,
         node_colors=node_colors,
         colors=colors
         )


class Test_Vertical_Line_Intersect (unittest.TestCase):

    def test_case(self):
        seg = Line(start=50 + 200j, end=50 + 700j)
        horiz1 = Line(start=0 + 200j, end=70 + 200j)
        horiz2 = Line(start=0 + 300j, end=70 + 300j)
        horiz3 = Line(start=0 + 700j, end=70 + 700j)
        intersections1 = seg.intersect(horiz1, tol=0.01)
        intersections2 = seg.intersect(horiz2, tol=0.01)
        intersections3 = seg.intersect(horiz3, tol=0.01)
        points1 = [seg.point(t1) for (t1, t2) in intersections1]
        points2 = [seg.point(t1) for (t1, t2) in intersections2]
        points3 = [seg.point(t1) for (t1, t2) in intersections3]
        #todisplay = sPath(seg, horiz1, horiz2, horiz3)
        #filename = self.id()[9:]
        #visualize_path(todisplay, filename, pathnodes=False, addpoints=[points1, points2, points3])
        self.assertEqual(len(points1), 1, msg="vertical line should intersect with its start point")
        self.assertEqual(len(points2), 1, msg="vertical line should intersect with at 1 point with middle line")
        self.assertEqual(len(points3), 1, msg="vertical line should intersect with its end point")
        self.assertAlmostEqual(points1[0], seg.start, msg="intersect point is not start point")
        self.assertAlmostEqual(points2[0], 50 + 300j, msg="intersect point with middle line is incorrect")
        self.assertAlmostEqual(points3[0], seg.end, msg="intersect point is not end point")


class Test_Curve_intersect():
    """Parent class of all curve test cases classes."""

    # Name used to generate the SVG: arc, quadractic_bezier or cubic_bezier
    name = ''

    # Tolerance under which a |derivative| will be considered zero
    ABSOLUTE_TOLERANCE = 0.01

    # List containing nominal cases
    paths = []

    # List containing edge cases
    edge_paths = []

    # list containing the results
    results = {}

    # specific ym list
    ym_list = None

    def generateSVG(self, path, i, addpoints, colors):
        """Create an SVG file containing test case number i."""
        filename = f'{self.name}_{i}'
        visualize_path(path, filename, pathnodes=False, addpoints=addpoints, colors=colors)
        self.pick_up_results(filename, addpoints, path[0])
        return filename

    def pick_up_results(self, filename, points, seg):
        """ this function picksup results which will be checked on in the test
        useful for non regression tests but visual checking is necessary during debug"""
        pointsS, pointsE, pointsM = points
        newPointsS = [p for p in pointsS if not isclose(p, seg.start, atol=self.ABSOLUTE_TOLERANCE)]
        newPointsE = [p for p in pointsE if not isclose(p, seg.end, atol=self.ABSOLUTE_TOLERANCE)]
        newPoints = [newPointsS, newPointsE, pointsM]
        countsS = len(pointsS)
        countsE = len(pointsE)
        countsM = len(pointsM)
        res = [countsS, countsE, countsM, newPoints]
        #print(f'"{filename}": {res},')

    def compute_lines(self, seg, ym):
        """ returns the 3 lines to intersect with seg"""
        ys = seg.start.imag
        ye = seg.end.imag
        if not ym:
            ym = (abs(ys - ye) / 2) + min(ys, ye)
        ps1 = complex(0, ys)
        ps2 = complex(max(p.real for p in [seg.start, seg.end]) + 500, ys)
        lineS = Line(start=ps1, end=ps2)
        pe1 = complex(0, ye)
        pe2 = complex(max(p.real for p in [seg.start, seg.end]) + 500, ye)
        lineE = Line(start=pe1, end=pe2)
        pm1 = complex(0, ym)
        pm2 = complex(max(p.real for p in [seg.start, seg.end]) + 500, ym)
        lineM = Line(start=pm1, end=pm2)
        return lineS, lineE, lineM

    def test_cases(self):
        """Check the results for all test cases."""
        # Looping over curves
        for i, path in enumerate(self.paths + self.edge_paths, start=1):
            seg = path[0]
            if self.ym_list:
                ym = self.ym_list[i - 1]
            else:
                ym = None
            lineS, lineE, lineM = self.compute_lines(seg, ym)
            intersectionsS = seg.intersect(lineS, tol=self.ABSOLUTE_TOLERANCE)
            intersectionsE = seg.intersect(lineE, tol=self.ABSOLUTE_TOLERANCE)
            intersectionsM = seg.intersect(lineM, tol=self.ABSOLUTE_TOLERANCE)
            pointsS = [seg.point(t1) for (t1, t2) in intersectionsS]
            pointsE = [seg.point(t1) for (t1, t2) in intersectionsE]
            pointsM = [seg.point(t1) for (t1, t2) in intersectionsM]
            todisplay = sPath(seg, lineS, lineE, lineM)
            colors = ('black', 'green', 'blue', 'black')
            filename = self.generateSVG(todisplay, i, addpoints=[pointsS, pointsE, pointsM], colors=colors)
            # for non regression test compare here the results with the values which were taken as output of the test
            with self.subTest(i=i, msg=filename):
                lenS, lenE, lenM, [respS, respE, respM] = self.results.get(filename)
                self.assertEqual(len(pointsS), lenS, msg="number of intersect with line at start(green) is incorrect")
                self.assertEqual(len(pointsE), lenE, msg='number of intersect with line at end (blue) is incorrect')
                self.assertEqual(len(pointsM), lenM, msg="number of intersect with line in the middle (black) is incorrect")
                # check start and end are is in the intersections
                startIsIn = any(isclose(point, seg.start, atol=self.ABSOLUTE_TOLERANCE) for point in pointsS)
                self.assertTrue(startIsIn, msg="line at start should intersect with segment")
                endIsIn = any(isclose(point, seg.end, atol=self.ABSOLUTE_TOLERANCE) for point in pointsE)
                self.assertTrue(endIsIn, msg='line at end should intersect with segment')
                # check that other points of intersections mentioned in results are correct
                for p in respS:
                    isin = any(isclose(point, p, atol=self.ABSOLUTE_TOLERANCE) for point in pointsS)
                    self.assertTrue(isin, msg=" other intersect with line at start(green) is missing")
                for p in respE:
                    isin = any(isclose(point, p, atol=self.ABSOLUTE_TOLERANCE) for point in pointsE)
                    self.assertTrue(isin, msg=" other intersect with line at end(blue) is missing")
                for p in respM:
                    isin = any(isclose(point, p, atol=self.ABSOLUTE_TOLERANCE) for point in pointsM)
                    self.assertTrue(isin, msg=" other intersect with line in the middle is missing")


class Test_Intersect_Arc (Test_Curve_intersect, unittest.TestCase):
    """Test case class for arcs."""

    @ classmethod
    def setUpClass(cls):
        """Initialize test cases for Arc class."""
        # configure_pathPoint_class()
        cls.name = 'arc'

        # Initialize  cases for Arc class
        for flag in ['0 0', '0 1', '1 0', '1 1']:
            cls.paths.extend([
                parse_path(f"M  50, 200 a 100, 50   0 {flag} 250,50"),
                parse_path(f"M 400, 100 a 100, 50  30 {flag} 250,50"),
                parse_path(f"M 400, 300 a 100, 50  45 {flag} 22,244"),
                parse_path(f"M 750, 200 a 100, 50 135 {flag} 22,244"),
                parse_path(f"M 750, 300 a 100, 50 220 {flag} 22,244"),
                parse_path(f"M 950, 100 a 100, 50 310 {flag} 22,244"),
            ])

        # Initialize edge cases for Arc class
        cls.edge_paths.extend([
            parse_path("M 400  200 a100 50, 0,  0, 0, 0 250"),
            parse_path("M 50 200 a100 50, 0,  0, 1, 0 250"),
            parse_path("M 400 200 a100 50, 90, 0, 0, 0 250"),
        ])
        # visual control of results show that most results are wrong if the arc has a rotation
        cls.results.update({
            "i1_arc_0_0_0": [1, 2, 1, [[], [(49.99999999999997 + 249.99999999999997j)], [(40.3708798216374 + 225j)]]],
            "i2_arc_30_0_0": [1, 2, 1, [[], [(437.1153753760908 + 149.9999999999999j)], []]],
            "i3_arc_45_0_0": [2, 1, 1, [[], [(275.59999999999997 + 300j)], [(295.26495776991305 + 421.999999999)]]],
            "i4_arc_135_0_0": [1, 2, 1, [[], [(603.5999999999999 + 443.9999999999994j)], [(632.0992629966765 + 321.9999999999999j)]]],
            "i5_arc_220_0_0": [2, 1, 1, [[(611.0555150486008 + 300j)], [], [(631.783782028338 + 421.9999999999999j)]]],
            "i6_arc_310_0_0": [1, 2, 1, [[], [(819.4282706707221 + 344.0000000000003j)], [(844.2419419089059 + 222.00000000000006j)]]],
            "i7_arc_0_0_1": [2, 1, 1, [[(300 + 200j)], [], [(309.6291201783626 + 225j)]]],
            "i8_arc_30_0_1": [2, 1, 1, [[(612.8846246239091 + 99.99999999999994j)], [], [(635.2102190245106 + 124.99999999999997j)]]],


            "i9_arc_45_0_1": [1, 2, 1, [[], [(546.3999999999999 + 544.0000000000002j)], [(526.735042230087 + 421.99999999999994j)]]],
            "i10_arc_135_0_1": [2, 1, 1, [[(918.4000000000001 + 200.00000000000034j)], [], [(889.9007370033237 + 321.99999999999994j), ]]],
            "i11_arc_220_0_1": [1, 2, 1, [[], [(910.9444849513991 + 544.0000000000001j)], [(890.2162179716619 + 421.99999999999994j), (631.783782028338 + 421.9999999999999j)]]],
            "i12_arc_310_0_1": [2, 1, 1, [[(1102.5717293292782 + 100.00000000000063j)], [(819.4282706707221 + 344.0000000000003j)], [(1077.758058091094 + 222.00000000000006j)]]],
            "i13_arc_0_1_0": [1, 2, 1, [[], [(49.99999999999997 + 249.99999999999997j)], [(40.3708798216374 + 225j)]]],
            "i14_arc_30_1_0": [1, 2, 1, [[], [(437.11537351971816 + 150.00000000000006j)], [(414.789780047303 + 124.99999999999993j)]]],
            "i15_arc_45_1_0": [2, 1, 1, [[(275.59999999999997 + 300j)], [], [(295.26495776991305 + 421.9999999999999j)]]],
            "i16_arc_135_1_0": [1, 2, 1, [[], [(603.5999999999999 + 443.9999999999994j)], [(632.0992629966765 + 321.9999999999999j)]]],
            "i17_arc_220_1_0": [2, 1, 1, [[(611.0555150486008 + 300j)], [], [(631.783782028338 + 421.9999999999999j)]]],
            "i18_arc_310_1_0": [1, 2, 1, [[], [(819.4282706707221 + 344.0000000000003j)], [(844.2419419089059 + 222.00000000000006j)]]],
            "i19_arc_0_1_1": [2, 1, 1, [[(300 + 200j)], [], [(309.6291201783626 + 225j)]]],
            "i20_arc_30_1_1": [2, 1, 1, [[(612.8846264802819 + 99.99999999999996j)], [], [(635.2102199526969 + 124.99999999999991j)]]],
            "i21_arc_45_1_1": [1, 2, 1, [[], [(546.3999999999999 + 544.0000000000002j)], [(526.735042230087 + 421.99999999999994j)]]],
            "i22_arc_135_1_1": [2, 1, 1, [[(918.4000000000001 + 200.00000000000034j)], [], [(889.9007370033237 + 321.99999999999994j)]]],
            "i23_arc_220_1_1": [1, 2, 1, [[], [(910.9444849513991 + 544.0000000000001j)], [(890.2162179716619 + 421.99999999999994j)]]],
            "i24_arc_310_1_1": [2, 1, 1, [[(1102.5717293292782 + 100.00000000000063j)], [], [(1077.758058091094 + 222.00000000000006j)]]],
            "i25_arc_0_0_0": [1, 1, 1, [[], [], [(150 + 325j)]]],
            "i26_arc_0_0_1": [1, 1, 1, [[], [], [(300 + 325j)]]],
            "i27_arc_90_0_0": [1, 1, 1, [[], [], [(337.5 + 325j)]]],
        })

    def generateSVG(self, path, i, addpoints, colors):
        """
        Override generateSVG for Arc class.

        The name of those files should contain additional information:
        rotation, large arc and sweep.
        """
        seg = path[0]
        filename = 'i%d_%s_%d_%d_%d' % (i, self.name, seg.rotation, seg.large_arc, seg.sweep)
        visualize_path(
            path,
            filename,
            pathnodes=False,
            addpoints=addpoints,
            colors=colors
        )
        # the following line is to generate the results dictionnary
        self.pick_up_results(filename, addpoints, seg)
        return filename

    def assertions(self, filename, pointsS, pointsE, pointsM):
        # these assertions are only for non regression tests. Presently the values in results are incorrect
        # all arcs with orientation give incorrect results
        res1, res2, res3 = self.results.get(filename)
        self.assertEqual(len(pointsS), len(res1))
        self.assertEqual(len(pointsE), len(res2))
        self.assertEqual(len(pointsM), len(res3))
        for index, p in enumerate(pointsS):
            self.assertAlmostEqual(p, res1[index])
        for index, p in enumerate(pointsE):
            self.assertAlmostEqual(p, res2[index])
        for index, p in enumerate(pointsM):
            self.assertAlmostEqual(p, res3[index])


class Test_Intersect_QuadraticBezier (Test_Curve_intersect, unittest.TestCase):
    """Test case class for quadratic bezier curves."""

    @ classmethod
    def setUpClass(cls):
        """Initialize test cases for Quadratic Bezier class."""
        cls.name = 'quadractic_bezier'

        # Initialize nominal cases for Quadratic Bezier class
        cls.paths = [
            parse_path("m  25, 47 q 17,  -7  6,  20"),
            parse_path("m  48, 68 q 14  -57 26,   0"),
            parse_path("m  85, 59 q 14,  57 26,   0"),
            parse_path("m 128, 64 q 14,  17 25, -18"),
            parse_path("m 128, 64 q 14,  17 25, -18"),
            parse_path("m  81, 29 q 14, -17 25,  18"),
        ]

        # Initialize edge cases for Quadratic Bezier class
        cls.edge_paths = [
            parse_path("m 80, 20 q 10,10 -10, 20"),
            # parse_path("m 80, 20 q  0, 0   0,  0"),  # Single point
            # parse_path("m 80, 20 q 10, 0  20,  0"),  # Straight line
            parse_path("m 80, 20 q 10,10   0, 20"),
        ]
        # ym_list for specific yms
        cls.ym_list = [None, 50, 80, None, None, None, None, None]

        # visual control all correct
        cls.results.update({
            "quadractic_bezier_1": [2, 1, 1, [[(34.25259515570934 + 46.99999999999999j)], [], [(34.425982139631 + 57j)]]],
            "quadractic_bezier_2": [2, 2, 2, [[(74 + 68j)], [(48 + 68j)], [(69.2064901963537 + 50j), (53.42508875101472 + 49.99999999999999j)]]],
            "quadractic_bezier_3": [2, 2, 2, [[(111 + 59j)], [(85 + 59j)], [(105.03728034118508 + 80j), (91.69956176407807 + 80j)]]],
            "quadractic_bezier_4": [2, 1, 1, [[(145.0251479289941 + 64j)], [], [(149.7705100077186 + 55j)]]],
            "quadractic_bezier_5": [2, 1, 1, [[(145.0251479289941 + 64j)], [], [(149.7705100077186 + 55j)]]],
            "quadractic_bezier_6": [2, 1, 1, [[(98.02514792899407 + 29j)], [], [(102.77051000771863 + 38j)]]],
            "quadractic_bezier_7": [1, 1, 1, [[], [], [(82.5 + 30j)]]],
            "quadractic_bezier_8": [1, 1, 1, [[], [], [(85 + 30j)]]],
        })


class Test_Intersect_CubicBezier (Test_Curve_intersect, unittest.TestCase):
    """Test case class for quadratic bezier curves."""

    @ classmethod
    def setUpClass(cls):
        """Initialize test cases for Cubic Bezier class."""
        cls.name = 'cubic_bezier'

        # Initialize nominal cases for Cubic Bezier class
        cls.paths = [
            parse_path("m  10, 110 c  46,-16  21,  22  5,  24"),
            parse_path("m  70, 134 c  15,  7  51, -18 29, -32"),
            parse_path("m 139, 100 c  26,-15  18,  33  8,  44"),
            parse_path("m  35, 151 c  15, -7  51,  18 29,  32"),
            parse_path("m 133, 154 c -46,-16 -21,  22 -5,  24"),
        ]

        # Initialize edge cases for Cubic Bezier class
        cls.edge_paths = [
            parse_path("m 139, 100 c 26, 5 18, 33  8, 44"),
            parse_path("m  35, 151 c 15, 0 51, 18 29, 32"),
            # parse_path("m  80,  20 c 10, 0 20,  0 30,  0"), this is a line
        ]
        # visual control cubic_bezier_1 and cubic_bezier_5 and cubic_bezier_7 are missing one cutpoint at the end
        cls.results.update({
            "cubic_bezier_1": [2, 1, 1, [[(36.003610559732806 + 110j)], [], [(32.078141538112355 + 122.00000000000001j)]]],
            "cubic_bezier_2": [2, 1, 1, [[(83.69506360982835 + 134j)], [], [(104.62766766084235 + 118j)]]],
            "cubic_bezier_3": [2, 1, 1, [[(153.52855800284362 + 100j)], [], [(155.59972699011112 + 121.99999999999997j)]]],
            "cubic_bezier_4": [2, 1, 1, [[(48.69506360982835 + 151j)], [], [(69.62766766084235 + 167j)]]],
            "cubic_bezier_5": [2, 1, 1, [[(106.9963894402672 + 154j)], [], [(110.92185846188764 + 166j)]]],
            "cubic_bezier_6": [1, 1, 1, [[], [], [(156.4095736229877 + 122j)]]],
            "cubic_bezier_7": [1, 1, 1, [[], [], [(68.65030968920104 + 167j)]]],
        })


if __name__ == '__main__':
    unittest.main()
