"""(Experimental) replacement for import/export functionality.

This module contains the `Document` class, a container for a DOM-style 
document (e.g. svg, html, xml, etc.) designed to replace and improve 
upon the IO functionality of svgpathtools (i.e. the svg2paths and 
disvg/wsvg functions). 

An Historic Note:
     The functionality in this module is meant to replace and improve 
     upon the IO functionality previously provided by the the 
     `svg2paths` and `disvg`/`wsvg` functions. 

Example:
    Typical usage looks something like the following.

        >> from svgpathtools import *
        >> doc = Document('my_file.html')
        >> for p in doc.paths:
        >>     foo(p)  # do stuff using svgpathtools functionality
        >> foo2(doc.tree)  # do stuff using ElementTree's functionality
        >> doc.display()  # display doc in OS's default application
        >> doc.save('my_new_file.html')

Todo: (please see contributor guidelines in CONTRIBUTING.md)
    * Finish "NotImplemented" methods.
    * Find some clever (and easy to implement) way to create a thorough 
    set of unittests.
    * Finish Documentation for each method (approximately following the 
    Google Python Style Guide, see [1]_ for some nice examples).  
    For nice style examples, see [1]_.

Some thoughts on this module's direction:
    * The `Document` class should ONLY grab path elements that are 
    inside an SVG.  
    * To handle transforms... there should be a "get_transform" 
    function and also a "flatten_transforms" tool that removes any 
    present transform attributes from all SVG-Path elements in the 
    document (applying the transformations before to the svgpathtools 
    Path objects).
    Note: This ability to "flatten" will ignore CSS files (and any 
    relevant files that are not parsed into the tree automatically by 
    ElementTree)... that is unless you have any bright ideas on this.  
    I really know very little about DOM-style documents.
"""

# External dependencies
from __future__ import division, absolute_import, print_function
import os
import xml.etree.cElementTree as etree

# Internal dependencies
from .parser import parse_path
from .svg2paths import (ellipse2pathd, line2pathd, polyline2pathd,
                        polygon2pathd, rect2pathd)
from .misctools import open_in_browser


CONVERSIONS = {'circle': ellipse2pathd,
               'ellipse': ellipse2pathd,
               'line': line2pathd,
               'polyline': polyline2pathd,
               'polygon': polygon2pathd,
               'rect': rect2pathd}


class Document:
    def __init__(self, filename=None, conversions=CONVERSIONS):
        """A container for a DOM-style document.

        The `Document` class is meant to be used to parse, create, save, 
        and modify DOM-style documents.  Given the `filename` of a 
        DOM-style document, it parses the document into an ElementTree 
        object, extracts all SVG-Path and Path-like (see `conversions` 
        below) objects into a list of svgpathtools Path objects."""

        # remember location of original svg file
        if filename is not None and os.path.dirname(filename) == '':
            self.original_filename = os.path.join(os.getcwd(), filename)
        else:
            self.original_filename = filename

        # parse svg to ElementTree object
        self.tree = etree.parse(filename)
        self.root = self.tree.getroot()

        # get URI namespace (only necessary in OS X?)
        root_tag = self.tree.getroot().tag
        if root_tag[0] == "{":
            self._prefix = root_tag[:root_tag.find('}') + 1]
        else:
            self._prefix = ''
        # etree.register_namespace('', prefix)

        self.paths = self._get_paths(conversions)

    def get_elements_by_tag(self, tag):
        """Returns a generator of all elements with the give tag. 

        Note: for more advanced XML-related functionality, use the 
        `tree` attribute (an ElementTree object).
        """
        return self.tree.iter(tag=self._prefix + tag)

    def _get_paths(self, conversions):
        paths = []

        # Get d-strings for SVG-Path elements
        paths += [el.attrib for el in self.get_elements_by_tag('path')]
        d_strings = [el['d'] for el in paths]
        attribute_dictionary_list = paths

        # Convert path-like elements to d-strings and attribute dicts
        if conversions:
            for tag, fcn in conversions.items():
                attributes = [el.attrib for el in self.get_elements_by_tag(tag)]
                d_strings += [fcn(d) for d in attributes]

        path_list = [parse_path(d) for d in d_strings]
        return path_list

    def get_svg_attributes(self):
        """To help with backwards compatibility."""
        return self.get_elements_by_tag('svg')[0].attrib

    def get_path_attributes(self):
        """To help with backwards compatibility."""
        return [p.tree_element.attrib for p in self.paths]

    def add(self, path, attribs={}, parent=None):
        """Add a new path to the SVG."""
        if parent is None:
            parent = self.tree.getroot()
        # just get root
        # then add new path
        # then record element_tree object in path
        raise NotImplementedError

    def add_group(self, group_attribs={}, parent=None):
        """Add an empty group element to the SVG."""
        if parent is None:
            parent = self.tree.getroot()
        raise NotImplementedError

    def update_tree(self):
        """Rewrite d-string's for each path in the `tree` attribute."""
        raise NotImplementedError

    def save(self, filename, update=True):
        """Write to svg to a file."""
        if update:
            self.update_tree()

        with open(filename, 'w') as output_svg:
            output_svg.write(etree.tostring(self.tree.getroot()))

    def display(self, filename=None, update=True):
        """Displays/opens the doc using the OS's default application."""
        if update:
            self.update_tree()

        if filename is None:
            raise NotImplementedError

        # write to a (by default temporary) file
        with open(filename, 'w') as output_svg:
            output_svg.write(etree.tostring(self.tree.getroot()))

        open_in_browser(filename)
