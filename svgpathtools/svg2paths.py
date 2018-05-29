"""A submodule of tools for creating path objects from SVG files.
The main tool being the svg2paths() function."""

# External dependencies
from __future__ import division, absolute_import, print_function
import os
import xml.etree.cElementTree as etree

# Internal dependencies
from .parser import parse_path


COORD_PAIR_TMPLT = re.compile(
    r'([\+-]?\d*[\.\d]\d*[eE][\+-]?\d+|[\+-]?\d*[\.\d]\d*)' +
    r'(?:\s*,\s*|\s+|(?=-))' +
    r'([\+-]?\d*[\.\d]\d*[eE][\+-]?\d+|[\+-]?\d*[\.\d]\d*)'
)


def ellipse2pathd(ellipse):
    """converts the parameters from an ellipse or a circle to a string 
    for a Path object d-attribute"""

    cx = ellipse.get('cx', None)
    cy = ellipse.get('cy', None)
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


def polyline2pathd(polyline_d, is_polygon=False):
    """converts the string from a polyline points-attribute to a string 
    for a Path object d-attribute"""
    points = COORD_PAIR_TMPLT.findall(polyline_d)
    closed = (float(points[0][0]) == float(points[-1][0]) and
              float(points[0][1]) == float(points[-1][1]))

    # The `parse_path` call ignores redundant 'z' (closure) commands
    # e.g. `parse_path('M0 0L100 100Z') == parse_path('M0 0L100 100L0 0Z')`
    # This check ensures that an n-point polygon is converted to an n-Line path.
    if is_polygon and closed:
        points.append(points[0])

    d = 'M' + 'L'.join('{0} {1}'.format(x,y) for x,y in points)
    if is_polygon or closed:
        d += 'z'
    return d


def polygon2pathd(polyline_d):
    """converts the string from a polygon points-attribute to a string 
    for a Path object d-attribute.
    Note:  For a polygon made from n points, the resulting path will be
    composed of n lines (even if some of these lines have length zero).
    """
    return polyline2pathd(polyline_d, True)


def rect2pathd(rect):
    """Converts an SVG-rect element to a Path d-string.

    The rectangle will start at the (x,y) coordinate specified by the 
    rectangle object and proceed counter-clockwise."""
    x0, y0 = float(rect.get('x', 0)), float(rect.get('y', 0))
    w, h = float(rect["width"]), float(rect["height"])
    x1, y1 = x0 + w, y0
    x2, y2 = x0 + w, y0 + h
    x3, y3 = x0, y0 + h

    d = ("M{} {} L {} {} L {} {} L {} {} z"
         "".format(x0, y0, x1, y1, x2, y2, x3, y3))
    return d


def line2pathd(l):
    return 'M' + l['x1'] + ' ' + l['y1'] + 'L' + l['x2'] + ' ' + l['y2']


CONVERSIONS = {'circle': ellipse2pathd,
               'ellipse': ellipse2pathd,
               'line': line2pathd,
               'polyline': polyline2pathd,
               'polygon': polygon2pathd,
               'rect': rect2pathd}


def svg2paths(svg_file_location, return_svg_attributes=False,
              conversions=CONVERSIONS, return_tree=False):
    """Converts SVG to list of Path objects and attribute dictionaries.

    Converts an SVG file into a list of Path objects and a list of
    dictionaries containing their attributes.  This currently supports
    SVG Path, Line, Polyline, Polygon, Circle, and Ellipse elements.

    Args:
        svg_file_location (string): the location of the svg file
        return_svg_attributes (bool): Set to True and a dictionary of
            svg-attributes will be extracted and returned.  See also 
            the `svg2paths2()` function.
        convert_circles_to_paths: Set to False to exclude SVG-Circle
            elements (converted to Paths).  By default circles are 
            included as paths of two `Arc` objects.
        convert_ellipses_to_paths (bool): Set to False to exclude 
            SVG-Ellipse elements (converted to Paths).  By default 
            ellipses are included as paths of two `Arc` objects.
        convert_lines_to_paths (bool): Set to False to exclude SVG-Line 
            elements (converted to Paths)
        convert_polylines_to_paths (bool): Set to False to exclude 
            SVG-Polyline elements (converted to Paths)
        convert_polygons_to_paths (bool): Set to False to exclude 
            SVG-Polygon elements (converted to Paths)
        convert_rectangles_to_paths (bool): Set to False to exclude 
            SVG-Rect elements (converted to Paths).

    Returns:
        list: The list of Path objects.
        list: The list of corresponding path attribute dictionaries.
        dict (optional): A dictionary of svg-attributes (see `
            svg2paths2()`).
    """
    if os.path.dirname(svg_file_location) == '':
        svg_file_location = os.path.join(getcwd(), svg_file_location)

    tree = etree.parse(svg_file_location)

    # get URI namespace
    root_tag = tree.getroot().tag
    if root_tag[0] == "{":
        prefix = root_tag[:root_tag.find('}') + 1]
    else:
        prefix = ''
    # etree.register_namespace('', prefix)

    def getElementsByTagName(tag):
        return tree.iter(tag=prefix+tag)

    # Get d-strings for Path elements
    paths = [el.attrib for el in getElementsByTagName('path')]
    d_strings = [el['d'] for el in paths]
    attribute_dictionary_list = paths

    # Get d-strings for Path-like elements (using `conversions` dict)
    for tag, fcn in conversions.items():
        attributes = [el.attrib for el in getElementsByTagName(tag)]
        d_strings += [fcn(d) for d in attributes]

    path_list = [parse_path(d) for d in d_strings]
    if return_tree:  # svg2paths3 default behavior
        return path_list, tree

    elif return_svg_attributes:  # svg2paths2 default behavior
        svg_attributes = getElementsByTagName('svg')[0].attrib
        return path_list, attribute_dictionary_list, svg_attributes

    else:  # svg2paths default behavior
        return path_list, attribute_dictionary_list


def svg2paths2(svg_file_location, return_svg_attributes=True,
               conversions=CONVERSIONS, return_tree=False):
    """Convenience function; identical to svg2paths() except that
    return_svg_attributes=True by default.  See svg2paths() docstring 
    for more info."""
    return svg2paths(svg_file_location=svg_file_location,
                     return_svg_attributes=return_svg_attributes,
                     conversions=conversions, return_tree=return_tree)


def svg2paths3(svg_file_location, return_svg_attributes=True,
               conversions=CONVERSIONS, return_tree=True):
    """Convenience function; identical to svg2paths() except that
    return_tree=True.  See svg2paths() docstring for more info."""
    return svg2paths(svg_file_location=svg_file_location,
                     return_svg_attributes=return_svg_attributes,
                     conversions=conversions, return_tree=return_tree)
