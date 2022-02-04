"""This submodule contains tools for creating path objects from SVG files.
The main tool being the svg2paths() function."""

# External dependencies
from __future__ import division, absolute_import, print_function
from xml.dom.minidom import parse
from os import path as os_path, getcwd
import re
from sys import version_info


# Internal dependencies
from .parser import parse_path


# f-strings are available in Python 3.6+
F_STRINGS_AVAILABLE = version_info < (3, 6)

COORD_PAIR_TMPLT = re.compile(
    r'([\+-]?\d*[\.\d]\d*[eE][\+-]?\d+|[\+-]?\d*[\.\d]\d*)' +
    r'(?:\s*,\s*|\s+|(?=-))' +
    r'([\+-]?\d*[\.\d]\d*[eE][\+-]?\d+|[\+-]?\d*[\.\d]\d*)'
)


def path2pathd(path):
    return path.get('d')


def ellipse2pathd(ellipse):
    """converts the parameters from an ellipse or a circle to a string for a 
    Path object d-attribute"""

    cx = ellipse.get('cx', 0)
    cy = ellipse.get('cy', 0)
    rx = ellipse.get('rx', None)
    ry = ellipse.get('ry', None)
    r = ellipse.get('r', None)

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
    """converts the string from a polyline points-attribute to a string for a
    Path object d-attribute"""
    if isinstance(polyline, str):
        points = polyline
    else:
        points = COORD_PAIR_TMPLT.findall(polyline.get('points', ''))

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
    """converts the string from a polygon points-attribute to a string 
    for a Path object d-attribute.
    Note:  For a polygon made from n points, the resulting path will be
    composed of n lines (even if some of these lines have length zero).
    """
    return polyline2pathd(polyline, True)


def rect2pathd(rect):
    """Converts an SVG-rect element to a Path d-string.
    
    The rectangle will start at the (x,y) coordinate specified by the 
    rectangle object and proceed counter-clockwise."""
    x, y = float(rect.get('x', 0)), float(rect.get('y', 0))
    w, h = float(rect.get('width', 0)), float(rect.get('height', 0))
    if 'rx' in rect or 'ry' in rect:

        # if only one, rx or ry, is present, use that value for both
        # https://developer.mozilla.org/en-US/docs/Web/SVG/Element/rect
        rx = rect.get('rx', None)
        ry = rect.get('ry', None)
        if rx is None:
            rx = ry or 0.
        if ry is None:
            ry = rx or 0.
        rx, ry = float(rx), float(ry)

        d = "M {} {} ".format(x + rx, y)  # right of p0
        d += "L {} {} ".format(x + w - rx, y)  # go to p1
        d += "A {} {} 0 0 1 {} {} ".format(rx, ry, x+w, y+ry)  # arc for p1
        d += "L {} {} ".format(x+w, y+h-ry)  # above p2
        d += "A {} {} 0 0 1 {} {} ".format(rx, ry, x+w-rx, y+h)  # arc for p2
        d += "L {} {} ".format(x+rx, y+h)  # right of p3
        d += "A {} {} 0 0 1 {} {} ".format(rx, ry, x, y+h-ry)  # arc for p3
        d += "L {} {} ".format(x, y+ry)  # below p0
        d += "A {} {} 0 0 1 {} {} z".format(rx, ry, x+rx, y)  # arc for p0
        return d

    x0, y0 = x, y
    x1, y1 = x + w, y
    x2, y2 = x + w, y + h
    x3, y3 = x, y + h

    d = ("M{} {} L {} {} L {} {} L {} {} z"
         "".format(x0, y0, x1, y1, x2, y2, x3, y3))
        
    return d


if F_STRINGS_AVAILABLE:
    def line2pathd(node_dict):
        return f"M{node_dict['x1']} {node_dict['y1']}L{node_dict['x2']} {node_dict['y2']}"
else:
    def line2pathd(node_dict):
        return "M{} {}L{} {}".format(node_dict['x1'], node_dict['y1'], node_dict['x2'], node_dict['y2'])


parser_dict = {
    'path': path2pathd,
    'circle': ellipse2pathd,
    'ellipse': ellipse2pathd,
    'line': line2pathd,
    'polyline': polyline2pathd,
    'polygon': polygon2pathd,
    'rect': rect2pathd,
}


def svg2paths(
    svg_file_location,
    return_svg_attributes=False,
    convert_circles_to_paths=True,
    convert_ellipses_to_paths=True,
    convert_lines_to_paths=True,
    convert_polylines_to_paths=True,
    convert_polygons_to_paths=True,
    convert_rectangles_to_paths=True,
):
    if os_path.dirname(svg_file_location) == "":
        svg_file_location = os_path.join(getcwd(), svg_file_location)

    doc = parse(svg_file_location)

    # todo: is this equivalent to `dict(element.attributes())`?
    def dom2dict(element):
        """Converts DOM elements to dictionaries of attributes."""
        keys = list(element.attributes.keys())
        values = [val.value for val in list(element.attributes.values())]
        d = dict(list(zip(keys, values)))
        # if not dict((k, v) for k, v in element.attributes.items()) == d:
        #     from IPython import embed; embed()  ### DEBUG
        return d

    include_dict = {
        'circle': convert_circles_to_paths,
        'ellipse': convert_ellipses_to_paths,
        'line': convert_lines_to_paths,
        'polyline': convert_polylines_to_paths,
        'polygon': convert_polygons_to_paths,
        'rect': convert_rectangles_to_paths,
    }
    included_elements = set(name for name, is_included in include_dict.items() if is_included)

    attribute_dictionaries = [
        (node.localName, dom2dict(node)) for node in doc.documentElement.childNodes
        if node.localName in included_elements
    ]
    d_strings = [parser_dict[element_type](nd) for element_type, nd in attribute_dictionaries]

    if return_svg_attributes:
        svg_attributes = dom2dict(doc.getElementsByTagName("svg")[0])
        doc.unlink()
        path_list = [parse_path(d) for d in d_strings]
        return path_list, attribute_dictionaries, svg_attributes
    else:
        doc.unlink()
        path_list = [parse_path(d) for d in d_strings]
        return path_list, attribute_dictionaries


def svg2paths2(svg_file_location,
               return_svg_attributes=True,
               convert_circles_to_paths=True,
               convert_ellipses_to_paths=True,
               convert_lines_to_paths=True,
               convert_polylines_to_paths=True,
               convert_polygons_to_paths=True,
               convert_rectangles_to_paths=True):
    """Convenience function; identical to svg2paths() except that
    return_svg_attributes=True by default.  See svg2paths() docstring for more
    info."""
    return svg2paths(svg_file_location=svg_file_location,
                     return_svg_attributes=return_svg_attributes,
                     convert_circles_to_paths=convert_circles_to_paths,
                     convert_ellipses_to_paths=convert_ellipses_to_paths,
                     convert_lines_to_paths=convert_lines_to_paths,
                     convert_polylines_to_paths=convert_polylines_to_paths,
                     convert_polygons_to_paths=convert_polygons_to_paths,
                     convert_rectangles_to_paths=convert_rectangles_to_paths)
