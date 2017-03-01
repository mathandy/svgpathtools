# External dependencies
from __future__ import division, absolute_import, print_function
from builtins import range
import unittest
from numpy import poly1d

# Internal dependencies
from svgpathtools import *


class TestPathTools(unittest.TestCase):

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
