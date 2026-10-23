import unittest
from os.path import exists, join
from tempfile import TemporaryDirectory
import xml.etree.ElementTree as ET

from svgpathtools import Line, Path, disvg, paths2Drawing, wsvg


class TestPaths2SVG(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = directory.name
        self.filename = join(self.directory, 'output.svg')
        self.path = Path(Line(0, 2.5 + 3.75j))

    def drawing(self, **kwargs):
        options = dict(filename=self.filename, stroke_widths=[0.1],
                       svgwrite_debug=True, openinbrowser=False)
        options.update(kwargs)
        return ET.fromstring(paths2Drawing(self.path, **options).tostring())

    def test_viewbox_dimensions_use_baseunit(self):
        cases = [({}, ('2.5', '3.75')),
                 ({'baseunit': 'px'}, ('2.5', '3.75')),
                 ({'baseunit': ''}, ('2.5', '3.75')),
                 ({'baseunit': None}, ('2.5', '3.75')),
                 ({'baseunit': 'cm'}, ('2.5cm', '3.75cm')),
                 ({'baseunit': 'mm'}, ('2.5mm', '3.75mm')),
                 ({'baseunit': 'in'}, ('2.5in', '3.75in')),
                 ({'baseunit': '%'}, ('2.5%', '3.75%'))]
        for viewbox in [(1, -2, 2.5, 3.75), '1 -2 2.5 3.75']:
            for options, expected in cases:
                with self.subTest(viewbox=viewbox, options=options):
                    root = self.drawing(viewbox=viewbox, **options)
                    self.assertEqual((root.get('width'), root.get('height')),
                                     expected)
                    self.assertEqual(root.get('viewBox'), '1 -2 2.5 3.75')
                    element = root.find('{http://www.w3.org/2000/svg}path')
                    self.assertEqual(element.get('d'), self.path.d())
        self.assertFalse(exists(self.filename))
        root = self.drawing(baseunit='cm', mindim=None, margin_size=0,
                            stroke_widths=[0])
        self.assertEqual((root.get('width'), root.get('height')),
                         ('2.5cm', '3.75cm'))

    def test_saved_viewbox_dimensions_use_baseunit(self):
        for writer in (disvg, wsvg):
            with self.subTest(writer=writer.__name__):
                filename = join(self.directory, writer.__name__ + '.svg')
                writer(self.path, filename=filename, openinbrowser=False,
                       stroke_widths=[0.1], viewbox=(1, -2, 2.5, 3.75),
                       baseunit='cm')
                root = ET.parse(filename).getroot()
                self.assertEqual((root.get('width'), root.get('height')),
                                 ('2.5cm', '3.75cm'))
                self.assertEqual(root.get('viewBox'), '1 -2 2.5 3.75')

    def test_explicit_dimensions_keep_their_units(self):
        for viewbox in (None, (1, -2, 2.5, 3.75)):
            for dimensions in (('7mm', '9mm'), (7, 9)):
                with self.subTest(viewbox=viewbox, dimensions=dimensions):
                    root = self.drawing(viewbox=viewbox,
                                        dimensions=dimensions, baseunit='cm')
                    self.assertEqual((root.get('width'), root.get('height')),
                                     tuple(map(str, dimensions)))

    def test_svg_attributes_override_inferred_dimensions(self):
        cases = [({}, ('2.5cm', '3.75cm')),
                 ({'width': '7in'}, ('7in', '3.75cm')),
                 ({'height': '9in'}, ('2.5cm', '9in')),
                 ({'width': '7in', 'height': '9in',
                   'viewBox': '0 0 7 9'}, ('7in', '9in'))]
        for attributes, expected in cases:
            with self.subTest(attributes=attributes):
                root = self.drawing(viewbox=(1, -2, 2.5, 3.75), baseunit='cm',
                                    svg_attributes=attributes)
                self.assertEqual((root.get('width'), root.get('height')),
                                 expected)
                if 'viewBox' in attributes:
                    self.assertEqual(root.get('viewBox'),
                                     attributes['viewBox'])
