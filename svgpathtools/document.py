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

A Big Problem:  
    Derivatives and other functions may be messed up by 
    transforms unless transforms are flattened (and not included in 
    css)
"""

# External dependencies
from __future__ import division, absolute_import, print_function
import os
import xml.etree.cElementTree as etree
import xml.etree.ElementTree.Element as Element
import xml.etree.ElementTree.SubElement as SubElement

# Internal dependencies
from .parser import parse_path
from .svg2paths import (ellipse2pathd, line2pathd, polyline2pathd,
                        polygon2pathd, rect2pathd)
from .misctools import open_in_browser
from .path import *

# THESE MUST BE WRAPPED TO OUPUT ElementTree.element objects
CONVERSIONS = {'circle': ellipse2pathd,
               'ellipse': ellipse2pathd,
               'line': line2pathd,
               'polyline': polyline2pathd,
               'polygon': polygon2pathd,
               'rect': rect2pathd}


def flatten_group_transforms(group):
    """Returns a 3x3 matrix which can transform points on a path from a group frame to the root frame"""
    if not isinstance(group, Element):
        raise TypeError('Must provide an xml.etree.Element object')




class Document:
    def __init__(self, filename, conversions=False, transform_paths=True):
        """(EXPERIMENTAL) A container for a DOM-style document.

        The `Document` class provides a simple interface to modify and analyze 
        the path elements in a DOM-style document.  The DOM-style document is 
        parsed into an ElementTree object (stored in the `tree` attribute) and
        all SVG-Path (and, optionally, Path-like) elements are extracted into a 
        list of svgpathtools Path objects. For more information on "Path-like"
        objects, see the below explanation of the `conversions` argument.
        
        Args:
            merge_transforms (object): 
            filename (str): The filename of the DOM-style object.
                conversions (bool or dict): If true, automatically converts 
                circle, ellipse, line, polyline, polygon, and rect elements 
                into path elements.  These changes are saved in the ElementTree 
                object.  For custom conversions, a dictionary can be passed in instead whose 
                keys are the element tags that are to be converted and whose values 
                are the corresponding conversion functions.  Conversion 
                functions should both take in and return an ElementTree.element
                object.
        """

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


    def get_elements_by_tag(self, tag):
        """Returns a generator of all elements with the given tag. 

        Note: for more advanced XML-related functionality, use the 
        `tree` attribute (an ElementTree object).
        """
        return self.tree.iter(tag=self._prefix + tag)

    def convert_pathlike_elements_to_paths(self, conversions=CONVERSIONS):
        raise NotImplementedError

    def get_svg_attributes(self):
        """To help with backwards compatibility."""
        return self.get_elements_by_tag('svg')[0].attrib

    def get_path_attributes(self):
        """To help with backwards compatibility."""
        return [p.tree_element.attrib for p in self.paths]

    def add_path(self, path, attribs={}, group=None):
        """Add a new path to the SVG."""

        # If we are not given a parent, assume that the path does not have a group
        if group is None:
            group = self.tree.getroot()

        # If we are given a list of strings (one or more), assume it represents a sequence of nested group names
        elif all(isinstance(elem, basestring) for elem in group):
            group = self.get_or_add_group(group)

        elif not isinstance(group, Element):
            raise TypeError('Must provide a list of strings or an xml.etree.Element object')

        # TODO: If the user passes in an xml.etree.Element object, should we check to make sure that it actually
        #       belongs to this Document object?

        if isinstance(path, Path):
            path_svg = path.d()
        elif is_path_segment(path):
            path_svg = Path(path).d()
        elif isinstance(path, basestring):
            # Assume this is a valid d-string TODO: Should we sanity check the input string?
            path_svg = path

        return SubElement(group, 'path', {'d': path_svg})

    def get_or_add_group(self, nested_names):
        """Get a group from the tree, or add a new one with the given name structure.

        *nested_names* is a list of strings which represent group names. Each group name will be nested inside of the
        previous group name.

        Returns the requested group. If the requested group did not exist, this function will create it, as well as all
        parent groups that it requires. All created groups will be left with blank attributes.

        """
        group = self.tree.getroot()
        # Drill down through the names until we find the desired group
        while nested_names:
            prev_group = group
            next_name = nested_names.pop(0)
            for elem in group.iter():
                if elem.get('id') == next_name:
                    group = elem

            if prev_group is group:
                # The group we're looking for does not exist, so let's create the group structure
                nested_names.insert(0, next_name)

                while nested_names:
                    next_name = nested_names.pop(0)
                    group = self.add_group({'id': next_name}, group)

                # Now nested_names will be empty, so the topmost while-loop will end

        return group

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
