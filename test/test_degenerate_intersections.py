import unittest

from svgpathtools import Arc, CubicBezier, Path, QuadraticBezier
from svgpathtools.bezier import bezier_intersections, boxes_intersect


class TestDegenerateIntersections(unittest.TestCase):
    def test_closed_bounding_boxes(self):
        cases = [
            ((0, 6, 0, 0), (2, 2, -4, 2), True),
            ((0, 6, 0, 0), (7, 7, -4, 2), False),
            ((0, 6, 0, 0), (2, 4, 0, 0), True),
            ((0, 6, 0, 0), (2, 4, 1, 1), False),
            ((0, 0, 0, 0), (-1, 1, -1, 1), True),
            ((0, 1, 0, 1), (1, 2, 1, 2), True),
            ((0, 1, 0, 1), (2, 3, 0, 1), False),
        ]
        for first, second, expected in cases:
            with self.subTest(first=first, second=second):
                self.assertEqual(boxes_intersect(first, second), expected)
                self.assertEqual(boxes_intersect(second, first), expected)

    def test_axis_aligned_beziers(self):
        horizontal = [QuadraticBezier(0, 3, 6), CubicBezier(0, 2, 4, 6),
                      QuadraticBezier(0, 0, 6), CubicBezier(0, 0, 0, 6)]
        vertical = [QuadraticBezier(2 - 4j, 2 - 1j, 2 + 2j),
                    CubicBezier(2 - 4j, 2 - 2j, 2, 2 + 2j),
                    CubicBezier(2 - 4j, 2 - 4j, 2 - 4j, 2 + 2j)]
        for first in horizontal:
            for second in vertical:
                for a, b in [(first, second), (second, first)]:
                    with self.subTest(a=a, b=b):
                        result = Path(a).intersect(Path(b))
                        self.assertEqual(len(result), 1)
                        (ta, _, _), (tb, _, _) = result[0]
                        self.assertAlmostEqual(a.point(ta), 2, delta=1e-10)
                        self.assertAlmostEqual(b.point(tb), 2, delta=1e-10)

    def test_flat_curve_against_curved_bezier(self):
        horizontal = CubicBezier(0, 2, 4, 6)
        # y(t) = 16*t**2 - 16*t + 3 has roots at 1/4 and 3/4.
        curved = CubicBezier(3j, 2 - 7j / 3, 4 - 7j / 3, 6 + 3j)
        result = sorted(horizontal.intersect(curved))
        self.assertEqual(len(result), 2)
        for (ta, tb), expected in zip(result, [1.5, 4.5]):
            self.assertAlmostEqual(horizontal.point(ta), expected, delta=1e-10)
            self.assertAlmostEqual(curved.point(tb), expected, delta=1e-10)

    def test_separated_flat_curves(self):
        horizontal = CubicBezier(0, 2, 4, 6)
        for other in [horizontal.translated(1j),
                      QuadraticBezier(7 - 4j, 7 - 1j, 7 + 2j)]:
            with self.subTest(other=other):
                self.assertEqual(horizontal.intersect(other), [])

    def test_endpoint_of_flat_curve(self):
        horizontal = CubicBezier(0, 2, 4, 6)
        vertical = CubicBezier(6 - 4j, 6 - 2j, 6, 6 + 2j)
        result = horizontal.intersect(vertical)
        self.assertEqual(len(result), 1)
        ta, tb = result[0]
        self.assertAlmostEqual(horizontal.point(ta), 6, delta=1e-10)
        self.assertAlmostEqual(vertical.point(tb), 6, delta=1e-10)

    def test_flat_bezier_against_arc(self):
        vertical = CubicBezier(-2j, -1j, 1j, 2j)
        arc = Arc(-1, 1 + 1j, 0, False, True, 1)
        for a, b in [(vertical, arc), (arc, vertical)]:
            with self.subTest(a=a, b=b):
                result = a.intersect(b)
                self.assertEqual(len(result), 1)
                ta, tb = result[0]
                self.assertAlmostEqual(a.point(ta), -1j, delta=1e-10)
                self.assertAlmostEqual(b.point(tb), -1j, delta=1e-10)

    def test_point_bezier_contacts(self):
        curve = CubicBezier(0, 2 + 2j, 4 - 2j, 6)
        for point, expected in [(3, True), (3 + 1j, False)]:
            node = CubicBezier(point, point, point, point)
            for a, b in [(node, curve), (curve, node)]:
                with self.subTest(a=a, b=b):
                    result = a.intersect(b)
                    self.assertEqual(len(result), int(expected))
                    for ta, tb in result:
                        self.assertAlmostEqual(a.point(ta), point)
                        self.assertAlmostEqual(b.point(tb), point)

    def test_point_bezier_on_flat_curve_or_arc(self):
        node = CubicBezier(-1j, -1j, -1j, -1j)
        for curve in [QuadraticBezier(-2j, 0, 2j),
                      Arc(-1, 1 + 1j, 0, False, True, 1)]:
            for a, b in [(node, curve), (curve, node)]:
                with self.subTest(a=a, b=b):
                    # Arc.intersect has its own polynomial special case;
                    # exercise the generic solver's Arc support directly.
                    result = (bezier_intersections(a, b, curve.length())
                              if isinstance(curve, Arc) else a.intersect(b))
                    self.assertEqual(len(result), 1)
                    ta, tb = result[0]
                    self.assertAlmostEqual(a.point(ta), -1j)
                    self.assertAlmostEqual(b.point(tb), -1j)

    def test_point_at_bezier_endpoint(self):
        curve = CubicBezier(35 + 151j, 50 + 151j, 86 + 169j, 64 + 183j)
        point = CubicBezier(curve.end, curve.end, curve.end, curve.end)
        self.assertEqual(curve.intersect(point), [(1., 0.)])
        self.assertEqual(point.intersect(curve), [(0., 1.)])

    def test_point_near_bezier_endpoint(self):
        curve = QuadraticBezier(0, 1, 2)
        point = QuadraticBezier(1e-9, 1e-9, 1e-9)
        for first, second in [(point, curve), (curve, point)]:
            with self.subTest(first=first):
                result = first.intersect(second, tol=1e-12)
                self.assertEqual(len(result), 1)
                t, u = result[0]
                self.assertLessEqual(
                    abs(first.point(t) - second.point(u)), 1e-12)

    def test_point_near_arc_endpoints(self):
        arc = Arc(-1, 1 + 2j, 37, False, True, 1)
        for curve in [arc, arc.reversed()]:
            for t in [1e-7, 1 - 1e-7]:
                z = curve.point(t)
                point = QuadraticBezier(z, z, z)
                for first, second in [(point, curve), (curve, point)]:
                    with self.subTest(first=first, t=t):
                        result = bezier_intersections(
                            first, second, curve.length(),
                            tol=1e-12, tol_deC=1e-12)
                        self.assertEqual(len(result), 1)
                        u, v = result[0]
                        self.assertLessEqual(
                            abs(first.point(u) - second.point(v)), 1e-12)

    def test_finite_collinear_contact_at_extrema(self):
        for direction in [1, 1j]:
            first = QuadraticBezier(*(p * direction + 3 + 4j
                                      for p in [0, 2, 0]))
            second = QuadraticBezier(*(p * direction + 3 + 4j
                                       for p in [2, 0, 2]))
            for a, b in [(first, second), (second, first)]:
                with self.subTest(a=a, b=b):
                    self.assertEqual(a.intersect(b), [(0.5, 0.5)])

    def test_finite_collinear_contact_at_endpoints(self):
        first = QuadraticBezier(0, 0.5, 1)
        second = CubicBezier(1, 2, 3, 4)
        self.assertEqual(first.intersect(second), [(1., 0.)])
        self.assertEqual(second.intersect(first), [(0., 1.)])


if __name__ == '__main__':
    unittest.main()
