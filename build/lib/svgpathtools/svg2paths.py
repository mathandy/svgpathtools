"""This submodule contains tools for creating path objects from SVG files.
The main tool being the svg2paths() function."""

# External dependencies
from __future__ import division, absolute_import, print_function
from xml.dom.minidom import parse
from os import path as os_path, getcwd
from shutil import copyfile

# Internal dependencies
from .parser import parse_path


def polyline2pathd(polyline_d):
    """converts the string from a polyline d-attribute to a string for a Path
    object d-attribute"""
    points = polyline_d.replace(', ', ',')
    points = points.replace(' ,', ',')
    points = points.split()

    if points[0] == points[-1]:
        closed = True
    else:
        closed = False

    d = 'M' + points.pop(0).replace(',', ' ')
    for p in points:
        d += 'L' + p.replace(',', ' ')
    if closed:
        d += 'z'
    return d


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

    def dom2dict(element):
        """Converts DOM elements to dictionaries of attributes."""
        keys = list(element.attributes.keys())
        values = [val.value for val in list(element.attributes.values())]
        return dict(list(zip(keys, values)))

    # Use minidom to extract path strings from input SVG
    paths = [dom2dict(el) for el in doc.getElementsByTagName('path')]
    d_strings = [el['d'] for el in paths]
    attribute_dictionary_list = paths
    # if pathless_svg:
    #     for el in doc.getElementsByTagName('path'):
    #         el.parentNode.removeChild(el)

    # Use minidom to extract polyline strings from input SVG, convert to
    # path strings, add to list
    if convert_polylines_to_paths:
        plins = [dom2dict(el) for el in doc.getElementsByTagName('polyline')]
        d_strings += [polyline2pathd(pl['points']) for pl in plins]
        attribute_dictionary_list += plins

    # Use minidom to extract polygon strings from input SVG, convert to
    # path strings, add to list
    if convert_polygons_to_paths:
        pgons = [dom2dict(el) for el in doc.getElementsByTagName('polygon')]
        d_strings += [polyline2pathd(pg['points']) + 'z' for pg in pgons]
        attribute_dictionary_list += pgons

    if convert_lines_to_paths:
        lines = [dom2dict(el) for el in doc.getElementsByTagName('line')]
        d_strings += [('M' + l['x1'] + ' ' + l['y1'] +
                       'L' + l['x2'] + ' ' + l['y2']) for l in lines]
        attribute_dictionary_list += lines

    # if pathless_svg:
    #     with open(pathless_svg, "wb") as f:
    #         doc.writexml(f)

    if return_svg_attributes:
        svg_attributes = dom2dict(doc.getElementsByTagName('svg')[0])
        doc.unlink()
        path_list = [parse_path(d) for d in d_strings]
        return path_list, attribute_dictionary_list, svg_attributes
    else:
        doc.unlink()
        path_list = [parse_path(d) for d in d_strings]
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
