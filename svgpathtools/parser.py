"""This submodule contains the path_parse() function used to convert SVG path
element d-strings into svgpathtools Path objects.
Note: This file was taken (nearly) as is from the svg.path module (v 2.0)."""

# External dependencies
from __future__ import division, absolute_import, print_function
import re
import numpy as np
import warnings

# Internal dependencies
from .path import Path, Line, QuadraticBezier, CubicBezier, Arc

# To maintain forward/backward compatibility
try:
    str = basestring
except NameError:
    pass


class PathTokens:
    """Path Tokens is the class for the general outline of how SVG Pathd objects
    are stored. Namely, a single non-'e' character and a collection of floating
    point numbers. While this is explicitly used for SVG pathd objects the method
    for serializing command data in this fashion is also useful as a standalone
    class."""

    def __init__(self, command_elements):
        self.command_elements = command_elements
        commands = ''
        for k in command_elements:
            commands += k
        self.COMMAND_RE = re.compile("([" + commands + "])")
        self.FLOAT_RE = re.compile("[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?")
        self.elements = None
        self.command = None
        self.last_command = None
        self.parser = None

    def _tokenize_path(self, pathdef):
        for x in self.COMMAND_RE.split(pathdef):
            if x in self.command_elements:
                yield x
            for token in self.FLOAT_RE.findall(x):
                yield token

    def get(self):
        """Gets the element from the stack."""
        return self.elements.pop()

    def pre_execute(self):
        """Called before any command element is executed."""
        pass

    def post_execute(self):
        """Called after any command element is executed."""
        pass

    def new_command(self):
        """Called when command element is switched."""
        pass

    def parse(self, pathdef):
        self.elements = list(self._tokenize_path(pathdef))
        # Reverse for easy use of .pop()
        self.elements.reverse()

        while self.elements:
            if self.elements[-1] in self.command_elements:
                self.last_command = self.command
                self.command = self.get()
                self.new_command()
            else:
                if self.command is None:
                    raise ValueError("Invalid command.")  # could be faulty implicit or unaccepted element.
            self.pre_execute()
            self.command_elements[self.command]()
            self.post_execute()


class SVGPathTokens(PathTokens):
    """Utilizes the general PathTokens class to parse SVG pathd strings.
    This class has been updated to account for SVG 2.0 version of the zZ command.
    Points are stored in complex numbers with the real being the x-value and
    the imaginary part being the y-value."""

    def __init__(self):
        PathTokens.__init__(self, {
            'M': self.move_to,
            'm': self.move_to,
            'L': self.line_to,
            'l': self.line_to,
            "H": self.h_to,
            "h": self.h_to,
            "V": self.v_to,
            "v": self.v_to,
            "C": self.cubic_to,
            "c": self.cubic_to,
            "S": self.smooth_cubic_to,
            "s": self.smooth_cubic_to,
            "Q": self.quad_to,
            "q": self.quad_to,
            "T": self.smooth_quad_to,
            "t": self.smooth_quad_to,
            "A": self.arc_to,
            "a": self.arc_to,
            "Z": self.close,
            "z": self.close
        })
        self.parser = None
        self.current_pos = None
        self.absolute = False
        self.start_pos = None
        self.last_cubic = None
        self.last_quad = None

    def svg_parse(self, parser, pathdef, current_pos=0j):
        # In the SVG specs, initial movetos are absolute, even if
        # specified as 'm'. This is the default behavior here as well.
        # But if you pass in a current_pos variable, the initial moveto
        # will be relative to that current_pos. This is useful.
        self.parser = parser
        self.current_pos = current_pos
        self.absolute = False
        self.start_pos = None
        self.last_cubic = None
        self.last_quad = None
        self.parser.start()
        self.parse(pathdef)
        self.parser.end()

    def get_pos(self):
        if self.command == 'Z':
            return self.start_pos  # After Z, all further expected values are also Z.
        coord0 = self.get()
        if coord0 == 'z' or coord0 == 'Z':
            self.command = 'Z'
            return self.start_pos
        coord1 = self.get()
        position = float(coord0) + float(coord1) * 1j
        if not self.absolute:
            return position + self.current_pos
        return position

    def move_to(self):
        # Moveto command.
        pos = self.get_pos()
        self.current_pos = pos

        # when M is called, reset start_pos
        # This behavior of Z is defined in svg spec:
        # http://www.w3.org/TR/SVG/paths.html#PathDataClosePathCommand
        self.start_pos = self.current_pos

        # Implicit moveto commands are treated as lineto commands.
        # So we set command to lineto here, in case there are
        # further implicit commands after this moveto.
        self.command = 'L'

    def line_to(self):
        pos = self.get_pos()
        self.parser.line(self.current_pos, pos)
        self.current_pos = pos

    def h_to(self):
        x = self.get()
        pos = float(x) + self.current_pos.imag * 1j
        if not self.absolute:
            pos += self.current_pos.real
        self.parser.line(self.current_pos, pos)
        self.current_pos = pos

    def v_to(self):
        y = self.get()
        pos = self.current_pos.real + float(y) * 1j
        if not self.absolute:
            pos += self.current_pos.imag * 1j
        self.parser.line(self.current_pos, pos)
        self.current_pos = pos

    def cubic_to(self):
        control1 = self.get_pos()
        control2 = self.get_pos()
        end = self.get_pos()
        self.parser.cubic(self.current_pos, control1, control2, end)
        self.current_pos = end
        self.last_cubic = control2

    def smooth_cubic_to(self):
        # Smooth curve. First control point is the "reflection" of
        # the second control point in the previous path.

        if self.last_command not in 'CScs':
            # If there is no previous command or if the previous command
            # was not an C, c, S or s, assume the first control point is
            # coincident with the current point.
            control1 = self.current_pos
        else:
            # The first control point is assumed to be the reflection of
            # the second control point on the previous command relative
            # to the current point.
            c_pos = self.current_pos
            control1 = c_pos + c_pos - self.last_cubic
        control2 = self.get_pos()
        end = self.get_pos()

        self.parser.cubic(self.current_pos, control1, control2, end)
        self.current_pos = end
        self.last_cubic = control2

    def quad_to(self):
        control = self.get_pos()
        end = self.get_pos()

        self.parser.quad(self.current_pos, control, end)
        self.current_pos = end
        self.last_quad = control

    def smooth_quad_to(self):
        # Smooth curve. Control point is the "reflection" of
        # the second control point in the previous path.

        if self.last_command not in 'QTqt':
            # If there is no previous command or if the previous command
            # was not an Q, q, T or t, assume the first control point is
            # coincident with the current point.
            control = self.current_pos
        else:
            # The control point is assumed to be the reflection of
            # the control point on the previous command relative
            # to the current point.
            c_pos = self.current_pos
            control = c_pos + c_pos - self.last_quad
        end = self.get_pos()

        self.parser.quad(self.current_pos, control, end)
        self.current_pos = end
        self.last_quad = control

    def arc_to(self):
        rx = self.get()
        ry = self.get()
        radius = float(rx) + float(ry) * 1j
        rotation = float(self.get())
        arc = float(self.get())
        sweep = float(self.get())
        end = self.get_pos()

        self.parser.arc(self.current_pos, radius, rotation, arc, sweep, end)
        self.current_pos = end

    def close(self):
        # Close path
        if not (self.current_pos == self.start_pos):
            self.parser.line(self.current_pos, self.start_pos)
        self.parser.closed()
        self.current_pos = self.start_pos
        self.command = None

    def new_command(self):
        self.absolute = self.command.isupper()

    def post_execute(self):
        if self.command == 'Z':  # Z might have been triggered inside commands.
            self.close()


class PathParser:
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

    def start(self):
        pass

    def end(self):
        pass

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
    tokens = SVGPathTokens()
    parser.set_path(tree_element)
    tokens.svg_parse(parser, pathdef, current_pos)
    return parser.segments


def _tokenize_transform(transform_str):
    if not transform_str:
        return
    transform_re = re.compile('(?u)(matrix|translate|scale|rotate|skewX|skewY)[\s\t\n]*\(([^)]+)\)')
    float_re = re.compile("[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?")
    for sub_element in transform_re.findall(transform_str):
        yield sub_element[0], list(map(float, float_re.findall(sub_element[1])))


def parse_svg_transform(transform_str, object):
    if not transform_str:
        return
    if not isinstance(transform_str, str):
        raise TypeError('Must provide a string to parse')

    object.identity()
    actions = list(_tokenize_transform(transform_str))
    for action in actions:
        if "matrix" in action[0]:
            object.matrix(action[1])
        elif "translate" in action[0]:
            object.translate(action[1])
        elif "scale" in action[0]:
            object.scale(action[1])
        elif "rotate" in action[0]:
            object.rotate(action[1])
        elif "skewX" in action[0]:
            object.skewX(action[1])
        elif "skewY" in action[0]:
            object.skewY(action[1])


class SvgMatrix:
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
