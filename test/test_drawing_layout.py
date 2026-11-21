import unittest

from svgpathtools import Line, paths2Drawing


class TestDrawingLayout(unittest.TestCase):
    def setUp(self):
        self.line = Line(0j, 1 + 1j)

    def test_explicit_viewbox_uses_default_styles(self):
        drawing = paths2Drawing(
            [self.line], nodes=[1j], viewbox=(0, 0, 10, 20)
        )
        self.assertEqual(drawing.attribs["viewBox"], "0 0 10 20")
        path, node = drawing.elements[1:]
        self.assertAlmostEqual(float(path.attribs["stroke-width"]), 0.02)
        self.assertAlmostEqual(float(node.attribs["r"]), 0.1)

    def test_explicit_dimensions_use_default_styles(self):
        drawing = paths2Drawing(
            [self.line], nodes=[1j], dimensions=("10mm", "20mm")
        )
        self.assertEqual(drawing.attribs["width"], "10mm")
        self.assertEqual(drawing.attribs["height"], "20mm")
        path, node = drawing.elements[1:]
        self.assertAlmostEqual(float(path.attribs["stroke-width"]), 0.02)
        self.assertAlmostEqual(float(node.attribs["r"]), 0.1)

    def test_explicit_styles_keep_their_values(self):
        drawing = paths2Drawing(
            [self.line],
            nodes=[1j],
            viewbox=(0, 0, 10, 20),
            attributes=[{"fill": "blue", "stroke": "green", "stroke-width": "2"}],
            node_radii=[3],
        )
        path, node = drawing.elements[1:]
        self.assertEqual(path.attribs["fill"], "blue")
        self.assertEqual(path.attribs["stroke"], "green")
        self.assertEqual(path.attribs["stroke-width"], "2")
        self.assertEqual(node.attribs["r"], 3)

    def test_explicit_stroke_width_is_independent_of_default_radius(self):
        drawing = paths2Drawing(
            [self.line], nodes=[1j], viewbox=(0, 0, 10, 20),
            stroke_widths=[0.75],
        )
        path, node = drawing.elements[1:]
        self.assertAlmostEqual(float(path.attribs["stroke-width"]), 0.75)
        self.assertAlmostEqual(float(node.attribs["r"]), 0.1)

    def test_automatic_layout_still_uses_geometry_scale(self):
        drawing = paths2Drawing([self.line], nodes=[1j])
        path, node = drawing.elements[1:]
        self.assertAlmostEqual(float(path.attribs["stroke-width"]), 0.001)
        self.assertAlmostEqual(float(node.attribs["r"]), 0.005)
