"""This submodule contains tools for creating path objects from SVG files.
The main tool being the svg2paths() function."""

# External dependencies
from __future__ import division, absolute_import, print_function
from xml.dom.minidom import parse
from os import path as os_path, getcwd
from shutil import copyfile
import numpy as np

# Internal dependencies
from .parser import parse_path
from .path import Path, bpoints2bezier


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
    # This check ensures that an n-point polygon is converted to an n-Line path.
    if reduntantly_closed:
        d += 'L' + points[0].replace(',', ' ')

    return d + 'z'


def svg2paths(svg_file_location,
              convert_lines_to_paths=True,
              convert_polylines_to_paths=True,
              convert_polygons_to_paths=True,
              return_svg_attributes=False):
    """
    Converts an SVG file into a list of Path objects and a list of
    dictionaries containing their attributes.  This currently supports
    SVG Path, Line, Polyline, and Polygon elements.
    :param svg_file_location: the location of the svg file
    :param convert_lines_to_paths: Set to False to disclude SVG-Line objects
    (converted to Paths)
    :param convert_polylines_to_paths: Set to False to disclude SVG-Polyline
    objects (converted to Paths)
    :param convert_polygons_to_paths: Set to False to disclude SVG-Polygon
    objects (converted to Paths)
    :param return_svg_attributes: Set to True and a dictionary of
    svg-attributes will be extracted and returned
    :return: list of Path objects, list of path attribute dictionaries, and
    (optionally) a dictionary of svg-attributes

    """
    if os_path.dirname(svg_file_location) == '':
        svg_file_location = os_path.join(getcwd(), svg_file_location)

    # if pathless_svg:
    #     copyfile(svg_file_location, pathless_svg)
    #     doc = parse(pathless_svg)
    # else:
    doc = parse(svg_file_location)

    # Parse a list of paths
    def dom2dict(element):
        """Converts DOM elements to dictionaries of attributes."""
        keys = list(element.attributes.keys())
        values = [val.value for val in list(element.attributes.values())]
        return dict(list(zip(keys, values)))

    def parseTrafo(trafoStr):
        """Returns six matrix elements for a matrix transformation for any valid SVG transformation string."""
        #print(trafoStr)
        valueStr = trafoStr.split('(')[1].split(')')[0]
        values = list(map(float, valueStr.split(',')))
        if 'translate' in trafoStr:
            x = values[0]
            y = values[1] if (len(values) > 1) else 0.
            return [1.,0.,0.,1.,x,y]
        elif 'scale' in trafoStr:
            x = values[0]
            y = values[1] if (len(values) > 1) else 0.
            return [x,0.,0.,y,0.,0.]
        elif 'rotate' in trafoStr:
            a = values[0]
            x = values[1] if (len(values) > 1) else 0.
            y = values[2] if (len(values) > 2) else 0.
            A = np.dot([cos(a),sin(a),-sin(a),cos(a),0.,0.,0.,0.,1.].reshape((3,3)),[1.,0.,0.,1.,-x,-y,0.,0.,1.].reshape((3,3)))
            A = list(np.dot([1.,0.,0.,1.,x,y,0.,0.,1.].reshape((3,3)),A).reshape((9,))[:6])
            return A
        elif 'skewX' in trafoStr:
            a = values[0]
            return [1.,0.,tan(a),1.,0.,0.]
        elif 'skewY' in trafoStr:
            a = values[0]
            return [1.,tan(a),0.,1.,0.,0.]
        else:
            while len(values) < 6:
               values += [0.]
            return values

    def parseNode(node):
        """Recursively iterate over nodes. Parse the groups individually to apply group transformations."""
        # Get everything in this tag
        #ret_list, attribute_dictionary_list = [parseNode(child) for child in node.childNodes]
        data = [parseNode(child) for child in node.childNodes]
        if len(data) == 0:
            ret_list = []
            attribute_dictionary_list = []
        else:
            # Flatten the lists
            ret_list = []
            attribute_dictionary_list = []
            for item in data:
                if type(item) == tuple:
                    if len(item[0]) > 0:
                        ret_list += item[0]
                        attribute_dictionary_list += item[1]
        
        if node.nodeName == 'g':
            # Group found
            # Analyse group properties
            group = dom2dict(node)
            if 'transform' in group.keys():
                trafo = group['transform']
                
                # Convert all transformations into a matrix operation
                A = parseTrafo(trafo)
                A = np.array([A[::2],A[1::2],[0.,0.,1.]])
                
                # Apply transformation to all elements of the paths
                xy = lambda z: np.array([z.real, z.imag, 1.])
                z = lambda xy: xy[0] + 1j*xy[1]
                
                ret_list = [Path(*[bpoints2bezier([z(np.dot(A,xy(pt)))
                    for pt in seg.bpoints()]) 
                        for seg in path])
                            for path in ret_list]
            return ret_list, attribute_dictionary_list
        elif node.nodeName == 'path':
            # Path found; parsing it
            path = dom2dict(node)
            d_string = path['d']
            return [parse_path(d_string)]+ret_list, [path]+attribute_dictionary_list
        elif convert_polylines_to_paths and node.nodeName == 'polyline':
            attrs = dom2dict(node)
            path = parse_path(polyline2pathd(node['points']))
            return [path]+ret_list, [attrs]+attribute_dictionary_list
        elif convert_polygons_to_paths and node.nodeName == 'polygon':
            attrs = dom2dict(node)
            path = parse_path(polygon2pathd(node['points']))
            return [path]+ret_list, [attrs]+attribute_dictionary_list
        elif convert_lines_to_paths and node.nodeName == 'line':
            line = dom2dict(node)
            d_string = ('M' + line['x1'] + ' ' + line['y1'] +
                        'L' + line['x2'] + ' ' + line['y2'])
            path = parse_path(d_string)
            return [path]+ret_list, [line]+attribute_dictionary_list
        else:
            return ret_list, attribute_dictionary_list

    path_list, attribute_dictionary_list = parseNode(doc)
    if return_svg_attributes:
        svg_attributes = dom2dict(doc.getElementsByTagName('svg')[0])
        doc.unlink()
        return path_list, attribute_dictionary_list, svg_attributes
    else:
        doc.unlink()
        return path_list, attribute_dictionary_list


def svg2paths2(svg_file_location,
               convert_lines_to_paths=True,
               convert_polylines_to_paths=True,
               convert_polygons_to_paths=True,
               return_svg_attributes=True):
    """Convenience function; identical to svg2paths() except that
    return_svg_attributes=True by default.  See svg2paths() docstring for more
    info."""
    return svg2paths(svg_file_location=svg_file_location,
                     convert_lines_to_paths=convert_lines_to_paths,
                     convert_polylines_to_paths=convert_polylines_to_paths,
                     convert_polygons_to_paths=convert_polygons_to_paths,
                     return_svg_attributes=return_svg_attributes)
