#!/usr/bin/env python
# coding:utf-8
"""
  SoftKnit21 framework

License
=======

 Copyright (C) <2019>  <Odile Troulet-Lamberttex@orange.fr>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or  any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.



Module purpose
==============

Tests for pathPoints

Implements
==========



Documentation
=============

Usage
=====


@author: Odile Troulet-Lambert
@copyright: Odile Troulet-Lambert
@license: GPL
"""
# from server.pathPoint import configure_pathPoint_class
from svgpathtools.path import Arc, QuadraticBezier, CubicBezier
from svgpathtools.parser import parse_path
from svgpathtools.paths2svg import wsvg
from pathlib import Path
import unittest


def visualize_path(path, label, addpoints=None, pathnodes=True):
    """Create an SVG file named label.svg containing path."""
    direct = './data/split_extremum_y_function'
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
        node_colors = ['red'] * len(startends) + ['blue'] * \
            len(arcs) + ['green'] * len(controls)
    else:
        nodes = []
        node_colors = []
    if addpoints and isinstance(addpoints[0], complex) and pathnodes:
        for point in addpoints:
            nodes.append(point)
        lenaddpoints = len(addpoints)
        node_colors += ['black'] * lenaddpoints
    elif addpoints and isinstance(addpoints[0], list) and not pathnodes:
        colors = ['red', 'blue', 'green', 'black']
        node_colors = []
        for i, lis in enumerate(addpoints):
            nodes += lis
            node_colors += [colors[i]] * len(lis)

    wsvg(path,
         filename=filename,
         nodes=nodes,
         node_colors=node_colors
         )


class Test_Split_extremum_y():
    """Parent class of all test cases classes."""

    # Name used to generate the SVG: arc, quadractic_bezier or cubic_bezier
    name = ''

    # Tolerance under which a |derivative| will be considered zero
    ABSOLUTE_TOLERANCE = 0.01

    # List containing nominal cases
    paths = []

    # List containing edge cases
    edge_paths = []

    def generateSVG(self, path, i, addpoints):
        """Create an SVG file containing test case number i."""
        visualize_path(path, f'{self.name}_{i}', addpoints=addpoints)

    def test_cases(self):
        """Check the results for all test cases."""
        # Looping over nominal test cases
        for i, path in enumerate(self.paths, start=1):
            seg = path[0]
            seg1, seg2 = seg.split_extremum_y()
            splitPoint = seg1.end
            self.generateSVG(path, i, [splitPoint])
            with self.subTest(i=i, msg=path):
                self.assertEqual(
                    splitPoint,
                    seg2.start,
                    msg="seg1.end should be equal to seg2.start"
                )
                self.assertAlmostEqual(
                    seg1.derivative(1).imag,
                    0,
                    delta=self.ABSOLUTE_TOLERANCE,
                    msg="The derivative should be 0 at the split point (left)"
                )
                self.assertAlmostEqual(
                    seg2.derivative(0).imag,
                    0,
                    delta=self.ABSOLUTE_TOLERANCE,
                    msg="The derivative should be 0 at the split point (right)"
                )

        # Looping over edge test cases
        for i, path in enumerate(self.edge_paths, start=len(self.paths) + 1):
            seg = path[0]
            seg1, seg2 = seg.split_extremum_y()
            self.generateSVG(path, i, None)
            with self.subTest(i=i, msg=path):
                self.assertTrue(
                    seg.__eq__(seg1),
                    msg="seg and seg1 should represent the same Curve"
                )
                self.assertIsNot(
                    seg,
                    seg1,
                    msg="seg1 should be a copy of seg"
                )
                self.assertIsNone(
                    seg2,
                    msg="seg2 should be None"
                )


class Test_SplitCurve_Arc (Test_Split_extremum_y, unittest.TestCase):
    """Test case class for arcs."""

    @classmethod
    def setUpClass(cls):
        """Initialize test cases for Arc class."""
#        configure_pathPoint_class()
        cls.name = 'arc'

        # Initialize nominal cases for Arc class
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
            parse_path("M 50 200 a100 50, 0,  0, 0, 0 250"),
            parse_path("M 50 200 a100 50, 0,  0, 1, 0 250"),
            parse_path("M 50 200 a100 50, 90, 0, 0, 0 250"),
        ])

    def generateSVG(self, path, i, addpoints):
        """
        Override generateSVG for Arc class.

        The name of those files should contain additional information:
        rotation, large arc and sweep.
        """
        seg = path[0]
        visualize_path(
            path,
            '%s_%d_%d_%d' % (self.name, seg.rotation, seg.large_arc, seg.sweep),
            addpoints=addpoints
        )


class Test_SplitCurve_QuadraticBezier (Test_Split_extremum_y, unittest.TestCase):
    """Test case class for quadratic bezier curves."""

    @classmethod
    def setUpClass(cls):
        """Initialize test cases for Quadratic Bezier class."""
        cls.name = 'quadractic_bezier'

        # Initialize nominal cases for Quadratic Bezier class
        cls.paths = [
            parse_path("m  25, 47 q 17,  -7  6,  20"),
            parse_path("m  48, 68 q 14  -17 26,   0"),
            parse_path("m  85, 59 q 14,  17 26,   0"),
            parse_path("m 128, 64 q 14,  17 25, -18"),
            parse_path("m 128, 64 q 14,  17 25, -18"),
            parse_path("m  81, 29 q 14, -17 25,  18"),
        ]

        # Initialize edge cases for Quadratic Bezier class
        cls.edge_paths = [
            parse_path("m 80, 20 q 10,10 -10, 20"),
            parse_path("m 80, 20 q  0, 0   0,  0"),  # Single point
            parse_path("m 80, 20 q 10, 0  20,  0"),  # Straight line
            parse_path("m 80, 20 q 10,10   0, 20"),
        ]


class Test_SplitCurve_CubicBezier (Test_Split_extremum_y, unittest.TestCase):
    """Test case class for quadratic bezier curves."""

    @classmethod
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
            parse_path("m  80,  20 c 10, 0 20,  0 30,  0"),
        ]


if __name__ == '__main__':
    unittest.main()
