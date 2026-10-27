import unittest

from svgpathtools import CubicBezier, Line, QuadraticBezier


class TestBezierLineEndpoints(unittest.TestCase):
    def test_reported_cubic_endpoints(self):
        examples = [
            (CubicBezier(10 + 110j, 56 + 94j, 31 + 132j, 15 + 134j),
             Line(134j, 515 + 134j)),
            (CubicBezier(133 + 154j, 87 + 138j, 112 + 176j, 128 + 178j),
             Line(178j, 633 + 178j)),
            (CubicBezier(35 + 151j, 50 + 151j, 86 + 169j, 64 + 183j),
             Line(183j, 564 + 183j)),
        ]
        for curve, line in examples:
            for reverse in [False, True]:
                a = curve.reversed() if reverse else curve
                endpoint = 0 if reverse else 1
                with self.subTest(curve=a, line=line):
                    result = a.intersect(line)
                    self.assertEqual(len(result), 1)
                    t, u = result[0]
                    self.assertEqual(t, endpoint)
                    self.assertAlmostEqual(a.point(t), line.point(u))

    def test_repeated_endpoint_with_interior_crossing(self):
        # y(t) = (t - 1/4)*(t - 1)**2 has a double endpoint root.
        curve = CubicBezier(-0.25j, 2 + 0.25j, 4, 6)
        line = Line(-1, 7)
        for a in [curve, curve.reversed()]:
            result = sorted(a.intersect(line))
            self.assertEqual(len(result), 2)
            expected = [0, 0.75] if a.start == curve.end else [0.25, 1]
            for (t, u), e in zip(result, expected):
                self.assertAlmostEqual(t, e)
                self.assertAlmostEqual(a.point(t), line.point(u))
            self.assertIn(0 if a.start == curve.end else 1,
                          [t for t, _ in result])

    def test_both_endpoints_and_interior(self):
        curve = CubicBezier(0, 2 - 1j / 6, 4 + 1j / 6, 6)
        line = Line(-1, 7)
        result = sorted(curve.intersect(line))
        self.assertEqual([t for t, _ in result], [0, 0.5, 1])
        for t, u in result:
            self.assertAlmostEqual(curve.point(t), line.point(u))

    def test_quadratic_tangent_at_endpoint(self):
        curve = QuadraticBezier(1j, 2, 4)
        line = Line(-1, 5)
        self.assertEqual(curve.intersect(line), [(1, 5 / 6)])
        self.assertEqual(curve.reversed().intersect(line), [(0, 5 / 6)])

    def test_nearby_curve_does_not_gain_an_endpoint(self):
        curve = CubicBezier(1j, 2 + 1j, 4 + 1j, 6 + 1j)
        self.assertEqual(curve.intersect(Line(-1, 7)), [])


if __name__ == '__main__':
    unittest.main()
