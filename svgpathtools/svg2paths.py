"""This submodule contains tools for creating path objects from SVG files.
The main tool being the svg2paths() function."""

# External dependencies
from __future__ import division, absolute_import, print_function
from xml.dom.minidom import parse
from os import path as os_path, getcwd
import numpy as np

# Internal dependencies
from .parser import parse_path
from .path import Path, bpoints2bezier


def ellipse2pathd(ellipse):
    """converts the parameters from an ellipse or a circle to a string for a 
    Path object d-attribute"""

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


def polyline2pathd(polyline_d):
    """converts the string from a polyline points-attribute to a string for a 
    Path object d-attribute"""
    points = polyline_d.replace(', ', ',')
    points = points.replace(' ,', ',')
    points = points.split()

    closed = points[0] == points[-1]

    d = 'M' + points.pop(0).replace(',', ' ')
    for p in points:
        d += 'L' + p.replace(',', ' ')
    if closed:
        d += 'z'
    return d


def polygon2pathd(polyline_d):
    """converts the string from a polygon points-attribute to a string for a 
    Path object d-attribute.
    Note:  For a polygon made from n points, the resulting path will be 
    composed of n lines (even if some of these lines have length zero)."""
    points = polyline_d.replace(', ', ',')
    points = points.replace(' ,', ',')
    points = points.split()

    reduntantly_closed = points[0] == points[-1]

    d = 'M' + points[0].replace(',', ' ')
    for p in points[1:]:
        d += 'L' + p.replace(',', ' ')

    # The `parse_path` call ignores redundant 'z' (closure) commands
    # e.g. `parse_path('M0 0L100 100Z') == parse_path('M0 0L100 100L0 0Z')`
    # This check ensures that an n-point polygon is converted to an n-Line
    # path.
    if reduntantly_closed:
        d += 'L' + points[0].replace(',', ' ')

    return d + 'z'


def rect2pathd(rect):
    """Converts an SVG-rect element to a Path d-string.
    
    The rectangle will start at the (x,y) coordinate specified by the rectangle 
    object and proceed counter-clockwise."""
    x0, y0 = float(rect.get('x', 0)), float(rect.get('y', 0))
    w, h = float(rect["width"]), float(rect["height"])
    x1, y1 = x0 + w, y0
    x2, y2 = x0 + w, y0 + h
    x3, y3 = x0, y0 + h

    d = ("M{} {} L {} {} L {} {} L {} {} z"
         "".format(x0, y0, x1, y1, x2, y2, x3, y3))
    return d


def svg2paths(svg_file_location,
              return_svg_attributes=False,
              convert_circles_to_paths=True,
              convert_ellipses_to_paths=True,
              convert_lines_to_paths=True,
              convert_polylines_to_paths=True,
              convert_polygons_to_paths=True,
              convert_rectangles_to_paths=True):
    """Converts an SVG into a list of Path objects and attribute dictionaries.

    Converts an SVG file into a list of Path objects and a list of
    dictionaries containing their attributes.  This currently supports
    SVG Path, Line, Polyline, Polygon, Circle, and Ellipse elements.

    Args:
        svg_file_location (string): the location of the svg file
        return_svg_attributes (bool): Set to True and a dictionary of
            svg-attributes will be extracted and returned.  See also the 
            `svg2paths2()` function.
        convert_circles_to_paths: Set to False to exclude SVG-Circle
            elements (converted to Paths).  By default circles are included as 
            paths of two `Arc` objects.
        convert_ellipses_to_paths (bool): Set to False to exclude SVG-Ellipse
            elements (converted to Paths).  By default ellipses are included as 
            paths of two `Arc` objects.
        convert_lines_to_paths (bool): Set to False to exclude SVG-Line elements
            (converted to Paths)
        convert_polylines_to_paths (bool): Set to False to exclude SVG-Polyline
            elements (converted to Paths)
        convert_polygons_to_paths (bool): Set to False to exclude SVG-Polygon
            elements (converted to Paths)
        convert_rectangles_to_paths (bool): Set to False to exclude SVG-Rect
            elements (converted to Paths).

    Returns: 
        list: The list of Path objects.
        list: The list of corresponding path attribute dictionaries.
        dict (optional): A dictionary of svg-attributes (see `svg2paths2()`).
    """
    if os_path.dirname(svg_file_location) == '':
        svg_file_location = os_path.join(getcwd(), svg_file_location)

    doc = parse(svg_file_location)

    def dom2dict(element):
        """Converts DOM elements to dictionaries of attributes."""
        keys = list(element.attributes.keys())
        values = [val.value for val in list(element.attributes.values())]
        return dict(list(zip(keys, values)))

    def parse_trafo(trafo_str):
        """Returns six matrix elements for a matrix transformation for any 
        valid SVG transformation string."""
        trafos = trafo_str.split(')')[:-1]
        trafo_matrix = np.array([1., 0., 0., 0., 1., 0., 0., 0., 1.]).reshape(
            (3, 3))  # Start with neutral matrix

        for trafo_sub_str in trafos:
            trafo_sub_str = trafo_sub_str.lstrip(', ')
            value_str = trafo_sub_str.split('(')[1]
            values = list(map(float, value_str.split(',')))
            if 'translate' in trafo_sub_str:
                x = values[0]
                y = values[1] if (len(values) > 1) else 0.
                trafo_matrix = np.dot(trafo_matrix, np.array(
                    [1., 0., x, 0., 1., y, 0., 0., 1.]).reshape((3, 3)))
            elif 'scale' in trafo_sub_str:
                x = values[0]
                y = values[1] if (len(values) > 1) else 0.
                trafo_matrix = np.dot(trafo_matrix,
                                      np.array([x, 0., 0., 0., y, 0., 0., 0.,
                                                1.]).reshape((3, 3)))
            elif 'rotate' in trafo_sub_str:
                a = values[0] * np.pi / 180.
                x = values[1] if (len(values) > 1) else 0.
                y = values[2] if (len(values) > 2) else 0.
                am = np.dot(np.array(
                    [np.cos(a), -np.sin(a), 0., np.sin(a), np.cos(a), 0., 0.,
                     0., 1.]).reshape((3, 3)),
                            np.array(
                                [1., 0., -x, 0., 1., -y, 0., 0., 1.]).reshape(
                                (3, 3)))
                am = np.dot(
                    np.array([1., 0., x, 0., 1., y, 0., 0., 1.]).reshape(
                        (3, 3)), am)
                trafo_matrix = np.dot(trafo_matrix, am)
            elif 'skewX' in trafo_sub_str:
                a = values[0] * np.pi / 180.
                trafo_matrix = np.dot(trafo_matrix,
                                      np.array(
                                          [1., np.tan(a), 0., 0., 1., 0., 0.,
                                           0., 1.]).reshape((3, 3)))
            elif 'skewY' in trafo_sub_str:
                a = values[0] * np.pi / 180.
                trafo_matrix = np.dot(trafo_matrix,
                                      np.array(
                                          [1., 0., 0., np.tan(a), 1., 0., 0.,
                                           0., 1.]).reshape((3, 3)))
            else:  # Assume matrix transformation
                while len(values) < 6:
                    values += [0.]
                trafo_matrix = np.dot(trafo_matrix,
                                      np.array([values[::2], values[1::2],
                                                [0., 0., 1.]]))

        trafo_list = list(trafo_matrix.reshape((9,))[:6])
        return trafo_list[::3] + trafo_list[1::3] + trafo_list[2::3]

    def parse_node(node):
        """Recursively iterate over nodes. Parse the groups individually to 
        apply group transformations."""
        # Get everything in this tag
        data = [parse_node(child) for child in node.childNodes]
        if len(data) == 0:
            ret_list = []
            attribute_dictionary_list_int = []
        else:
            # Flatten the lists
            ret_list = []
            attribute_dictionary_list_int = []
            for item in data:
                if type(item) == tuple:
                    if len(item[0]) > 0:
                        ret_list += item[0]
                        attribute_dictionary_list_int += item[1]

        if node.nodeName == 'g':
            # Group found
            # Analyse group properties
            group = dom2dict(node)
            if 'transform' in group.keys():
                trafo = group['transform']

                # Convert all transformations into a matrix operation
                am = parse_trafo(trafo)
                am = np.array([am[::2], am[1::2], [0., 0., 1.]])

                # Apply transformation to all elements of the paths
                def xy(p):
                    return np.array([p.real, p.imag, 1.])

                def z(coords):
                    return coords[0] + 1j * coords[1]

                ret_list = [Path(*[bpoints2bezier([z(np.dot(am, xy(pt)))
                                                   for pt in seg.bpoints()])
                                   for seg in path])
                            for path in ret_list]
            return ret_list, attribute_dictionary_list_int
        elif node.nodeName == 'path':
            # Path found; parsing it
            path = dom2dict(node)
            d_string = path['d']
            return [parse_path(d_string)] + ret_list, [
                path] + attribute_dictionary_list_int
        elif convert_polylines_to_paths and node.nodeName == 'polyline':
            attrs = dom2dict(node)
            path = parse_path(polyline2pathd(node['points']))
            return [path] + ret_list, [attrs] + attribute_dictionary_list_int
        elif convert_polygons_to_paths and node.nodeName == 'polygon':
            attrs = dom2dict(node)
            path = parse_path(polygon2pathd(attrs['points']))
            return [path] + ret_list, [attrs] + attribute_dictionary_list_int
        elif convert_lines_to_paths and node.nodeName == 'line':
            line = dom2dict(node)
            d_string = ('M' + line['x1'] + ' ' + line['y1'] +
                        'L' + line['x2'] + ' ' + line['y2'])
            path = parse_path(d_string)
            return [path] + ret_list, [line] + attribute_dictionary_list_int
        elif convert_ellipses_to_paths and node.nodeName == 'ellipse':
            attrs = dom2dict(node)
            path = parse_path(ellipse2pathd(attrs))
            return [path] + ret_list, [attrs] + attribute_dictionary_list_int
		elif convert_circles_to_paths and node.nodeName == 'circle':
			attrs = dom2dict(node)
            path = parse_path(ellipse2pathd(attrs))
            return [path] + ret_list, [attrs] + attribute_dictionary_list_int
		elif convert_rectangles_to_paths and node.nodeName == 'rect':
            attrs = dom2dict(node)
            path = parse_path(rect2pathd(attrs))
            return [path] + ret_list, [attrs] + attribute_dictionary_list_int
        else:
            return ret_list, attribute_dictionary_list_int

    path_list, attribute_dictionary_list = parse_node(doc)
    if return_svg_attributes:
        svg_attributes = dom2dict(doc.getElementsByTagName('svg')[0])
        doc.unlink()
        return path_list, attribute_dictionary_list, svg_attributes
    else:
        doc.unlink()
        return path_list, attribute_dictionary_list


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
