"""
(Experimental) replacement for import/export functionality SAX
"""

import re
from xml.etree.ElementTree import iterparse

# SVG STATIC VALUES
SVG_NAME_TAG = 'svg'
SVG_ATTR_VERSION = 'version'
SVG_VALUE_VERSION = '1.1'
SVG_ATTR_XMLNS = 'xmlns'
SVG_VALUE_XMLNS = 'http://www.w3.org/2000/svg'
SVG_ATTR_XMLNS_LINK = 'xmlns:xlink'
SVG_VALUE_XLINK = 'http://www.w3.org/1999/xlink'
SVG_ATTR_XMLNS_EV = 'xmlns:ev'
SVG_VALUE_XMLNS_EV = 'http://www.w3.org/2001/xml-events'

SVG_ATTR_WIDTH = 'width'
SVG_ATTR_HEIGHT = 'height'
SVG_ATTR_VIEWBOX = 'viewBox'
SVG_TAG_PATH = 'path'
SVG_TAG_GROUP = 'g'
SVG_TAG_RECT = 'rect'
SVG_TAG_CIRCLE = 'circle'
SVG_TAG_ELLIPSE = 'ellipse'
SVG_TAG_LINE = 'line'
SVG_TAG_POLYLINE = 'polyline'
SVG_TAG_POLYGON = 'polygon'
SVG_ATTR_DATA = 'd'
SVG_ATTR_FILL = 'fill'
SVG_ATTR_STROKE = 'stroke'
SVG_ATTR_STROKE_WIDTH = 'stroke-width'
SVG_ATTR_TRANSFORM = 'transform'
SVG_ATTR_STYLE = 'style'
SVG_ATTR_CENTER_X = 'cx'
SVG_ATTR_CENTER_Y = 'cy'
SVG_ATTR_RADIUS_X = 'rx'
SVG_ATTR_RADIUS_Y = 'ry'
SVG_ATTR_RADIUS = 'r'
SVG_ATTR_POINTS = 'points'
SVG_ATTR_X = 'x'
SVG_ATTR_Y = 'y'
SVG_ATTR_TAG = 'tag'
SVG_TRANSFORM_MATRIX = 'matrix'
SVG_TRANSFORM_TRANSLATE = 'translate'
SVG_TRANSFORM_SCALE = 'scale'
SVG_TRANSFORM_ROTATE = 'rotate'
SVG_TRANSFORM_SKEW_X = 'skewX'
SVG_TRANSFORM_SKEW_Y = 'skewY'
SVG_VALUE_NONE = 'none'

COORD_PAIR_TMPLT = re.compile(
    r'([\+-]?\d*[\.\d]\d*[eE][\+-]?\d+|[\+-]?\d*[\.\d]\d*)' +
    r'(?:\s*,\s*|\s+|(?=-))' +
    r'([\+-]?\d*[\.\d]\d*[eE][\+-]?\d+|[\+-]?\d*[\.\d]\d*)'
)


# Leaf node to pathd values.

def path2pathd(path):
    return path.get(SVG_ATTR_DATA, '')


def ellipse2pathd(ellipse):
    """converts the parameters from an ellipse or a circle to a string for a
    Path object d-attribute"""

    cx = ellipse.get(SVG_ATTR_CENTER_X, None)
    cy = ellipse.get(SVG_ATTR_CENTER_Y, None)
    rx = ellipse.get(SVG_ATTR_RADIUS_X, None)
    ry = ellipse.get(SVG_ATTR_RADIUS_X, None)
    r = ellipse.get(SVG_ATTR_RADIUS, None)

    if r is not None:
        rx = ry = float(r)
    else:
        rx = float(rx)
        ry = float(ry)

    cx = float(cx)
    cy = float(cy)

    d = ''
    d += 'M' + str(cx - rx) + ',' + str(cy)
    d += 'a' + str(rx) + ',' + str(ry) + ' 0 1,0 ' + str(2 * rx) + ',0'
    d += 'a' + str(rx) + ',' + str(ry) + ' 0 1,0 ' + str(-2 * rx) + ',0'

    return d


def polyline2pathd(polyline, is_polygon=False):
    """converts the string from a polyline parameters to a string for a
    Path object d-attribute"""
    polyline_d = polyline.get(SVG_ATTR_POINTS, None)
    if polyline_d is None:
        return ''
    points = COORD_PAIR_TMPLT.findall(polyline_d)
    closed = (float(points[0][0]) == float(points[-1][0]) and
              float(points[0][1]) == float(points[-1][1]))

    # The `parse_path` call ignores redundant 'z' (closure) commands
    # e.g. `parse_path('M0 0L100 100Z') == parse_path('M0 0L100 100L0 0Z')`
    # This check ensures that an n-point polygon is converted to an n-Line path.
    if is_polygon and closed:
        points.append(points[0])

    d = 'M' + 'L'.join('{0} {1}'.format(x, y) for x, y in points)
    if is_polygon or closed:
        d += 'z'
    return d


def polygon2pathd(polyline):
    """converts the string from a polygon parameters to a string
    for a Path object d-attribute.
    Note:  For a polygon made from n points, the resulting path will be
    composed of n lines (even if some of these lines have length zero).
    """
    return polyline2pathd(polyline, True)


def rect2pathd(rect):
    """Converts an SVG-rect element to a Path d-string.

    The rectangle will start at the (x,y) coordinate specified by the
    rectangle object and proceed counter-clockwise."""
    x0, y0 = float(rect.get(SVG_ATTR_X, 0)), float(rect.get(SVG_ATTR_Y, 0))
    w, h = float(rect.get(SVG_ATTR_WIDTH, 0)), float(rect.get(SVG_ATTR_HEIGHT, 0))
    x1, y1 = x0 + w, y0
    x2, y2 = x0 + w, y0 + h
    x3, y3 = x0, y0 + h

    d = ("M{} {} L {} {} L {} {} L {} {} z"
         "".format(x0, y0, x1, y1, x2, y2, x3, y3))
    return d


def line2pathd(l):
    return 'M' + l['x1'] + ' ' + l['y1'] + 'L' + l['x2'] + ' ' + l['y2']


# PathTokens class.
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


# Default Parser interface
class DefaultParser:
    """Default Parser gives an example of the needed thing without actually implementing anything"""

    def start(self):
        pass

    def end(self):
        pass

    def move(self, start_pos):
        pass

    def line(self, start_pos, end_pos):
        pass

    def quad(self, start_pos, control, end_pos):
        pass

    def cubic(self, start_pos, control1, control2, end_pos):
        pass

    def arc(self, start_pos, radius, rotation, arc, sweep, end_pos):
        pass

    def closed(self):
        pass


# SVG Path Tokens.
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
        self.parser.move(self.start_pos)

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


def parse_svg_path(parser, pathdef, current_pos=0j):
    """Parses the SVG path.
    The parser is datastructure agnostic requires a subclass of DefaultParser"""
    tokens = SVGPathTokens()
    tokens.svg_parse(parser, pathdef, current_pos)


class DefaultTransform:
    """Default Transform gives an example of the needed class without actually implementing anything"""

    def identity(self):
        pass

    def matrix(self, values):
        pass

    def translate(self, values):
        pass

    def scale(self, values):
        pass

    def rotate(self, values):
        pass

    def skewX(self, values):
        pass

    def skewY(self, values):
        pass


def _tokenize_transform(transform_str):
    """Generator to create transform parse elements.
    Will return tuples(command, list(values))"""
    if not transform_str:
        return
    transform_regex = '(?u)(' \
                      + SVG_TRANSFORM_MATRIX + '|' \
                      + SVG_TRANSFORM_TRANSLATE + '|' \
                      + SVG_TRANSFORM_SCALE + '|' \
                      + SVG_TRANSFORM_ROTATE + '|' \
                      + SVG_TRANSFORM_SKEW_X + '|' \
                      + SVG_TRANSFORM_SKEW_Y + \
                      ')[\s\t\n]*\(([^)]+)\)'
    transform_re = re.compile(transform_regex)
    float_re = re.compile("[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?")
    for sub_element in transform_re.findall(transform_str):
        yield sub_element[0], list(map(float, float_re.findall(sub_element[1])))


def parse_svg_transform(transform_str, obj):
    """Parses the svg transform tag.
    The object is datastructure agnostic. And requires a subclass of DefaultTransform"""
    if not transform_str:
        return
    if not isinstance(transform_str, str):
        raise TypeError('Must provide a string to parse')

    obj.identity()
    actions = list(_tokenize_transform(transform_str))
    for action in actions:
        if SVG_TRANSFORM_MATRIX in action[0]:
            obj.matrix(action[1])
        elif SVG_TRANSFORM_TRANSLATE in action[0]:
            obj.translate(action[1])
        elif SVG_TRANSFORM_SCALE in action[0]:
            obj.scale(action[1])
        elif SVG_TRANSFORM_ROTATE in action[0]:
            obj.rotate(action[1])
        elif SVG_TRANSFORM_SKEW_X in action[0]:
            obj.skewX(action[1])
        elif SVG_TRANSFORM_SKEW_Y in action[0]:
            obj.skewY(action[1])


def parse_svg_file(f):
    """Parses the SVG file.
    Style elements are split into their proper values.
    Transform elements are concatenated and unparsed.
    Leaf node elements are turned into pathd values."""

    stack = []
    values = {}
    for event, elem in iterparse(f, events=('start', 'end')):
        if event == 'start':
            stack.append(values)

            current_values = values
            values = {}
            values.update(current_values)  # copy of dictionary

            attrs = elem.attrib
            if SVG_ATTR_STYLE in attrs:
                for equate in attrs[SVG_ATTR_STYLE].split(";"):
                    equal_item = equate.split(":")
                    attrs[equal_item[0]] = equal_item[1]
            if SVG_ATTR_TRANSFORM in attrs:
                current_transform = values.get(SVG_ATTR_TRANSFORM, "")
                attrs[SVG_ATTR_TRANSFORM] = current_transform + attrs[SVG_ATTR_TRANSFORM]

            values.update(attrs)
            tag = elem.tag[28:]  # Removing namespace. http://www.w3.org/2000/svg:
            if SVG_NAME_TAG == tag:
                yield values
                continue
            elif SVG_TAG_GROUP == tag:
                continue
            elif SVG_TAG_PATH == tag:
                values[SVG_ATTR_DATA] = path2pathd(values)
            elif SVG_TAG_CIRCLE == tag:
                values[SVG_ATTR_DATA] = ellipse2pathd(values)
            elif SVG_TAG_ELLIPSE == tag:
                values[SVG_ATTR_DATA] = ellipse2pathd(values)
            elif SVG_TAG_LINE == tag:
                values[SVG_ATTR_DATA] = line2pathd(values)
            elif SVG_TAG_POLYLINE == tag:
                values[SVG_ATTR_DATA] = polyline2pathd(values)
            elif SVG_TAG_POLYGON == tag:
                values[SVG_ATTR_DATA] = polygon2pathd(values)
            elif SVG_TAG_RECT == tag:
                values[SVG_ATTR_DATA] = rect2pathd(values)
            else:
                continue
            values[SVG_ATTR_TAG] = tag
            yield values
        else:
            values = stack.pop()


# SVG Color Parsing

def color_rgb(r, g, b):
    return int(0xFF000000 |
               ((r & 255) << 16) |
               ((g & 255) << 8) |
               (b & 255))


# defining predefined colors permitted by svg: https://www.w3.org/TR/SVG11/types.html#ColorKeywords
svg_color_dict = {
    "aliceblue": color_rgb(240, 248, 255),
    "antiquewhite": color_rgb(250, 235, 215),
    "aqua": color_rgb(0, 255, 255),
    "aquamarine": color_rgb(127, 255, 212),
    "azure": color_rgb(240, 255, 255),
    "beige": color_rgb(245, 245, 220),
    "bisque": color_rgb(255, 228, 196),
    "black": color_rgb(0, 0, 0),
    "blanchedalmond": color_rgb(255, 235, 205),
    "blue": color_rgb(0, 0, 255),
    "blueviolet": color_rgb(138, 43, 226),
    "brown": color_rgb(165, 42, 42),
    "burlywood": color_rgb(222, 184, 135),
    "cadetblue": color_rgb(95, 158, 160),
    "chartreuse": color_rgb(127, 255, 0),
    "chocolate": color_rgb(210, 105, 30),
    "coral": color_rgb(255, 127, 80),
    "cornflowerblue": color_rgb(100, 149, 237),
    "cornsilk": color_rgb(255, 248, 220),
    "crimson": color_rgb(220, 20, 60),
    "cyan": color_rgb(0, 255, 255),
    "darkblue": color_rgb(0, 0, 139),
    "darkcyan": color_rgb(0, 139, 139),
    "darkgoldenrod": color_rgb(184, 134, 11),
    "darkgray": color_rgb(169, 169, 169),
    "darkgreen": color_rgb(0, 100, 0),
    "darkgrey": color_rgb(169, 169, 169),
    "darkkhaki": color_rgb(189, 183, 107),
    "darkmagenta": color_rgb(139, 0, 139),
    "darkolivegreen": color_rgb(85, 107, 47),
    "darkorange": color_rgb(255, 140, 0),
    "darkorchid": color_rgb(153, 50, 204),
    "darkred": color_rgb(139, 0, 0),
    "darksalmon": color_rgb(233, 150, 122),
    "darkseagreen": color_rgb(143, 188, 143),
    "darkslateblue": color_rgb(72, 61, 139),
    "darkslategray": color_rgb(47, 79, 79),
    "darkslategrey": color_rgb(47, 79, 79),
    "darkturquoise": color_rgb(0, 206, 209),
    "darkviolet": color_rgb(148, 0, 211),
    "deeppink": color_rgb(255, 20, 147),
    "deepskyblue": color_rgb(0, 191, 255),
    "dimgray": color_rgb(105, 105, 105),
    "dimgrey": color_rgb(105, 105, 105),
    "dodgerblue": color_rgb(30, 144, 255),
    "firebrick": color_rgb(178, 34, 34),
    "floralwhite": color_rgb(255, 250, 240),
    "forestgreen": color_rgb(34, 139, 34),
    "fuchsia": color_rgb(255, 0, 255),
    "gainsboro": color_rgb(220, 220, 220),
    "ghostwhite": color_rgb(248, 248, 255),
    "gold": color_rgb(255, 215, 0),
    "goldenrod": color_rgb(218, 165, 32),
    "gray": color_rgb(128, 128, 128),
    "grey": color_rgb(128, 128, 128),
    "green": color_rgb(0, 128, 0),
    "greenyellow": color_rgb(173, 255, 47),
    "honeydew": color_rgb(240, 255, 240),
    "hotpink": color_rgb(255, 105, 180),
    "indianred": color_rgb(205, 92, 92),
    "indigo": color_rgb(75, 0, 130),
    "ivory": color_rgb(255, 255, 240),
    "khaki": color_rgb(240, 230, 140),
    "lavender": color_rgb(230, 230, 250),
    "lavenderblush": color_rgb(255, 240, 245),
    "lawngreen": color_rgb(124, 252, 0),
    "lemonchiffon": color_rgb(255, 250, 205),
    "lightblue": color_rgb(173, 216, 230),
    "lightcoral": color_rgb(240, 128, 128),
    "lightcyan": color_rgb(224, 255, 255),
    "lightgoldenrodyellow": color_rgb(250, 250, 210),
    "lightgray": color_rgb(211, 211, 211),
    "lightgreen": color_rgb(144, 238, 144),
    "lightgrey": color_rgb(211, 211, 211),
    "lightpink": color_rgb(255, 182, 193),
    "lightsalmon": color_rgb(255, 160, 122),
    "lightseagreen": color_rgb(32, 178, 170),
    "lightskyblue": color_rgb(135, 206, 250),
    "lightslategray": color_rgb(119, 136, 153),
    "lightslategrey": color_rgb(119, 136, 153),
    "lightsteelblue": color_rgb(176, 196, 222),
    "lightyellow": color_rgb(255, 255, 224),
    "lime": color_rgb(0, 255, 0),
    "limegreen": color_rgb(50, 205, 50),
    "linen": color_rgb(250, 240, 230),
    "magenta": color_rgb(255, 0, 255),
    "maroon": color_rgb(128, 0, 0),
    "mediumaquamarine": color_rgb(102, 205, 170),
    "mediumblue": color_rgb(0, 0, 205),
    "mediumorchid": color_rgb(186, 85, 211),
    "mediumpurple": color_rgb(147, 112, 219),
    "mediumseagreen": color_rgb(60, 179, 113),
    "mediumslateblue": color_rgb(123, 104, 238),
    "mediumspringgreen": color_rgb(0, 250, 154),
    "mediumturquoise": color_rgb(72, 209, 204),
    "mediumvioletred": color_rgb(199, 21, 133),
    "midnightblue": color_rgb(25, 25, 112),
    "mintcream": color_rgb(245, 255, 250),
    "mistyrose": color_rgb(255, 228, 225),
    "moccasin": color_rgb(255, 228, 181),
    "navajowhite": color_rgb(255, 222, 173),
    "navy": color_rgb(0, 0, 128),
    "oldlace": color_rgb(253, 245, 230),
    "olive": color_rgb(128, 128, 0),
    "olivedrab": color_rgb(107, 142, 35),
    "orange": color_rgb(255, 165, 0),
    "orangered": color_rgb(255, 69, 0),
    "orchid": color_rgb(218, 112, 214),
    "palegoldenrod": color_rgb(238, 232, 170),
    "palegreen": color_rgb(152, 251, 152),
    "paleturquoise": color_rgb(175, 238, 238),
    "palevioletred": color_rgb(219, 112, 147),
    "papayawhip": color_rgb(255, 239, 213),
    "peachpuff": color_rgb(255, 218, 185),
    "peru": color_rgb(205, 133, 63),
    "pink": color_rgb(255, 192, 203),
    "plum": color_rgb(221, 160, 221),
    "powderblue": color_rgb(176, 224, 230),
    "purple": color_rgb(128, 0, 128),
    "red": color_rgb(255, 0, 0),
    "rosybrown": color_rgb(188, 143, 143),
    "royalblue": color_rgb(65, 105, 225),
    "saddlebrown": color_rgb(139, 69, 19),
    "salmon": color_rgb(250, 128, 114),
    "sandybrown": color_rgb(244, 164, 96),
    "seagreen": color_rgb(46, 139, 87),
    "seashell": color_rgb(255, 245, 238),
    "sienna": color_rgb(160, 82, 45),
    "silver": color_rgb(192, 192, 192),
    "skyblue": color_rgb(135, 206, 235),
    "slateblue": color_rgb(106, 90, 205),
    "slategray": color_rgb(112, 128, 144),
    "slategrey": color_rgb(112, 128, 144),
    "snow": color_rgb(255, 250, 250),
    "springgreen": color_rgb(0, 255, 127),
    "steelblue": color_rgb(70, 130, 180),
    "tan": color_rgb(210, 180, 140),
    "teal": color_rgb(0, 128, 128),
    "thistle": color_rgb(216, 191, 216),
    "tomato": color_rgb(255, 99, 71),
    "turquoise": color_rgb(64, 224, 208),
    "violet": color_rgb(238, 130, 238),
    "wheat": color_rgb(245, 222, 179),
    "white": color_rgb(255, 255, 255),
    "whitesmoke": color_rgb(245, 245, 245),
    "yellow": color_rgb(255, 255, 0),
    "yellowgreen": color_rgb(154, 205, 50)
}


def parse_svg_color_lookup(color_string):
    """Parse SVG Color by Keyword on dictionary lookup"""
    return svg_color_dict.get(color_string, 0xFF000000)


def parse_svg_color_hex(hex_string):
    """Parse SVG Color by Hex String"""
    h = hex_string.lstrip('#')
    size = len(h)
    if size == 8:
        return int(h[:8], 16)
    elif size == 6:
        return int(h[:6], 16)
    elif size == 4:
        return int(h[3] + h[3] + h[2] + h[2] + h[1] + h[1] + h[0] + h[0], 16)
    elif size == 3:
        return int(h[2] + h[2] + h[1] + h[1] + h[0] + h[0], 16)
    return 0xFF000000


def parse_svg_color_rgb(values):
    """Parse SVG Color, RGB value declarations """
    int_values = list(map(int, values))
    return color_rgb(int_values[0], int_values[1], int_values[2])


def parse_svg_color_rgbp(values):
    """Parse SVG color, RGB percent value declarations"""
    ratio = 255.0 / 100.0
    values = list(map(float, values))
    return color_rgb(int(values[0] * ratio), int(values[1] * ratio), int(values[2] * ratio))


def parse_svg_color(color_string):
    """Parse SVG color, will return a set value."""
    hex_re = re.compile(r'^#?([0-9A-Fa-f]{3,8})$')
    match = hex_re.match(color_string)
    if match:
        return parse_svg_color_hex(color_string)
    rgb_re = re.compile(r'rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)')
    match = rgb_re.match(color_string)
    if match:
        return parse_svg_color_rgb(match.groups())

    rgbp_re = re.compile(r'rgb\(\s*(\d+)%\s*,\s*(\d+)%\s*,\s*(\d+)%\s*\)')
    match = rgbp_re.match(color_string)
    if match:
        return parse_svg_color_rgbp(match.groups())
    return parse_svg_color_lookup(color_string)
