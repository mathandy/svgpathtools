"""This submodule contains the path_parse() function used to convert SVG path
element d-strings into svgpathtools Path objects.
Note: This file was taken (nearly) as is from the svg.path module (v 2.0)."""

# External dependencies
from __future__ import division, absolute_import, print_function
import numpy as np
import warnings

# Internal dependencies
from .path import Path, Line, QuadraticBezier, CubicBezier, Arc
from .svg_parser import DefaultParser, DefaultTransform, parse_svg_path, parse_svg_transform

# To maintain forward/backward compatibility
try:
    str = basestring
except NameError:
    pass


class PathParser(DefaultParser):
    """Path parser object for passing to the SVGPathTokens class to build svgpathtools/Path objects.
    This is intended to decouple the building of svgpathtool/Path objects from the rest of the parsing."""

    def __init__(self, tree_element=None):
        self.segments = None
        self.set_path(tree_element)

    def set_path(self, tree_element=None):
        if tree_element is None:
            self.segments = Path()
        else:
            self.segments = Path(tree_element)

    def line(self, start_pos, end_pos):
        self.segments.append(Line(start_pos, end_pos))

    def quad(self, start_pos, control, end_pos):
        self.segments.append(QuadraticBezier(start_pos, control, end_pos))

    def cubic(self, start_pos, control1, control2, end_pos):
        self.segments.append(CubicBezier(start_pos, control1, control2, end_pos))

    def arc(self, start_pos, radius, rotation, arc, sweep, end_pos):
        self.segments.append(Arc(start_pos, radius, rotation, arc, sweep, end_pos))

    def closed(self):
        self.segments.closed = True


def parse_path(pathdef, current_pos=0j, tree_element=None):
    parser = PathParser()
    parser.set_path(tree_element)
    parse_svg_path(parser, pathdef, current_pos)
    return parser.segments


class SvgMatrix(DefaultTransform):
    def __init__(self):
        self.transform = np.identity(3)

    @staticmethod
    def _check_num_parsed_values(values, allowed):
        if not any(num == len(values) for num in allowed):
            if len(allowed) > 1:
                warnings.warn('Expected one of the following number of values {0}, but found {1} values instead: {2}'
                              .format(allowed, len(values), values))
            elif allowed[0] != 1:
                warnings.warn('Expected {0} values, found {1}: {2}'.format(allowed[0], len(values), values))
            else:
                warnings.warn('Expected 1 value, found {0}: {1}'.format(len(values), values))
            return False
        return True

    def identity(self):
        self.transform = np.identity(3)

    def matrix(self, values):
        if not self._check_num_parsed_values(values, [6]):
            return

        element = np.identity(3)
        element[0:2, 0:3] = np.array([values[0:6:2], values[1:6:2]])
        self.transform = self.transform.dot(element)

    def translate(self, values):
        if not self._check_num_parsed_values(values, [1, 2]):
            return
        element = np.identity(3)
        element[0, 2] = values[0]
        if len(values) > 1:
            element[1, 2] = values[1]
        self.transform = self.transform.dot(element)

    def scale(self, values):
        if not self._check_num_parsed_values(values, [1, 2]):
            return
        element = np.identity(3)
        x_scale = values[0]
        y_scale = values[1] if (len(values) > 1) else x_scale
        element[0, 0] = x_scale
        element[1, 1] = y_scale
        self.transform = self.transform.dot(element)

    def rotate(self, values):
        if not self._check_num_parsed_values(values, [1, 3]):
            return
        angle = values[0] * np.pi / 180.0
        if len(values) == 3:
            offset = values[1:3]
        else:
            offset = (0, 0)
        tf_offset = np.identity(3)
        tf_offset[0:2, 2:3] = np.array([[offset[0]], [offset[1]]])
        tf_rotate = np.identity(3)
        tf_rotate[0:2, 0:2] = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        tf_offset_neg = np.identity(3)
        tf_offset_neg[0:2, 2:3] = np.array([[-offset[0]], [-offset[1]]])

        element = tf_offset.dot(tf_rotate).dot(tf_offset_neg)
        self.transform = self.transform.dot(element)

    def skewX(self, values):
        if not self._check_num_parsed_values(values, [1]):
            return
        element = np.identity(3)
        element[0, 1] = np.tan(values[0] * np.pi / 180.0)
        self.transform = self.transform.dot(element)

    def skewY(self, values):
        if not self._check_num_parsed_values(values, [1]):
            return
        element = np.identity(3)
        element[1, 0] = np.tan(values[0] * np.pi / 180.0)
        self.transform = self.transform.dot(element)


def parse_transform(transform_str):
    """Converts a valid SVG transformation string into a numpy 3x3 matrix.
        If the string is empty or null, this returns a 3x3 identity matrix"""
    matrix = SvgMatrix()
    parse_svg_transform(transform_str, matrix)
    transform = matrix.transform
    return transform
