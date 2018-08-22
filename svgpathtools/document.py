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
        >> results = doc.flatten_all_paths()
        >> for result in results:
        >>     path = result.path
        >>     # Do something with the transformed Path object.
        >>     element = result.element
        >>     # Inspect the raw SVG element. This gives access to the
        >>     # path's attributes
        >>     transform = result.transform
        >>     # Use the transform that was applied to the path.
        >> foo(doc.tree)  # do stuff using ElementTree's functionality
        >> doc.display()  # display doc in OS's default application
        >> doc.save('my_new_file.html')

A Big Problem:  
    Derivatives and other functions may be messed up by 
    transforms unless transforms are flattened (and not included in 
    css)
"""

# External dependencies
from __future__ import division, absolute_import, print_function
import os
import collections
import xml.etree.ElementTree as etree
from xml.etree.ElementTree import Element, SubElement, register_namespace
import warnings

# Internal dependencies
from .parser import parse_path
from .parser import parse_transform
from .svg_to_paths import (path2pathd, ellipse2pathd, line2pathd,
                           polyline2pathd, polygon2pathd, rect2pathd)
from .misctools import open_in_browser
from .path import *

# To maintain forward/backward compatibility
try:
    str = basestring
except NameError:
    pass

# Let xml.etree.ElementTree know about the SVG namespace
SVG_NAMESPACE = {'svg': 'http://www.w3.org/2000/svg'}
register_namespace('svg', 'http://www.w3.org/2000/svg')

# THESE MUST BE WRAPPED TO OUTPUT ElementTree.element objects
CONVERSIONS = {'path': path2pathd,
               'circle': ellipse2pathd,
               'ellipse': ellipse2pathd,
               'line': line2pathd,
               'polyline': polyline2pathd,
               'polygon': polygon2pathd,
               'rect': rect2pathd}

CONVERT_ONLY_PATHS = {'path': path2pathd}

SVG_GROUP_TAG = 'svg:g'


def flatten_all_paths(group, group_filter=lambda x: True,
                      path_filter=lambda x: True, path_conversions=CONVERSIONS,
                      group_search_xpath=SVG_GROUP_TAG):
    """Returns the paths inside a group (recursively), expressing the
    paths in the base coordinates.

    Note that if the group being passed in is nested inside some parent
    group(s), we cannot take the parent group(s) into account, because
    xml.etree.Element has no pointer to its parent. You should use
    Document.flatten_group(group) to flatten a specific nested group into
    the root coordinates.

    Args:
        group is an Element
        path_conversions (dict):
            A dictionary to convert from an SVG element to a path data
            string. Any element tags that are not included in this
            dictionary will be ignored (including the `path` tag). To
            only convert explicit path elements, pass in
            `path_conversions=CONVERT_ONLY_PATHS`.
    """
    if not isinstance(group, Element):
        raise TypeError('Must provide an xml.etree.Element object. '
                        'Instead you provided {0}'.format(type(group)))

    # Stop right away if the group_selector rejects this group
    if not group_filter(group):
        return []

    # To handle the transforms efficiently, we'll traverse the tree of
    # groups depth-first using a stack of tuples.
    # The first entry in the tuple is a group element and the second
    # entry is its transform. As we pop each entry in the stack, we
    # will add all its child group elements to the stack.
    StackElement = collections.namedtuple('StackElement',
                                          ['group', 'transform'])

    def new_stack_element(element, last_tf):
        return StackElement(element, last_tf.dot(
            parse_transform(element.get('transform'))))

    def get_relevant_children(parent, last_tf):
        children = []
        for elem in filter(group_filter,
                           parent.iterfind(group_search_xpath, SVG_NAMESPACE)):
            children.append(new_stack_element(elem, last_tf))
        return children

    stack = [new_stack_element(group, np.identity(3))]

    FlattenedPath = collections.namedtuple('FlattenedPath',
                                           ['path', 'element', 'transform'])
    paths = []

    while stack:
        top = stack.pop()

        # For each element type that we know how to convert into path
        # data, parse the element after confirming that the path_filter
        # accepts it.
        for key, converter in path_conversions.items():
            for path_elem in filter(path_filter, top.group.iterfind(
                    'svg:'+key, SVG_NAMESPACE)):
                path_tf = top.transform.dot(
                    parse_transform(path_elem.get('transform')))
                path = transform(parse_path(converter(path_elem)), path_tf)
                paths.append(FlattenedPath(path, path_elem, path_tf))

        stack.extend(get_relevant_children(top.group, top.transform))

    return paths


def flatten_group(group_to_flatten, root, recursive=True,
                  group_filter=lambda x: True, path_filter=lambda x: True,
                  path_conversions=CONVERSIONS,
                  group_search_xpath=SVG_GROUP_TAG):
    """Flatten all the paths in a specific group.

    The paths will be flattened into the 'root' frame. Note that root
    needs to be an ancestor of the group that is being flattened.
    Otherwise, no paths will be returned."""

    if not any(group_to_flatten is descendant for descendant in root.iter()):
        warnings.warn('The requested group_to_flatten is not a '
                      'descendant of root')
        # We will shortcut here, because it is impossible for any paths
        # to be returned anyhow.
        return []

    # We create a set of the unique IDs of each element that we wish to
    # flatten, if those elements are groups. Any groups outside of this
    # set will be skipped while we flatten the paths.
    desired_groups = set()
    if recursive:
        for group in group_to_flatten.iter():
            desired_groups.add(id(group))
    else:
        desired_groups.add(id(group_to_flatten))

    def desired_group_filter(x):
        return (id(x) in desired_groups) and group_filter(x)

    return flatten_all_paths(root, desired_group_filter, path_filter,
                             path_conversions, group_search_xpath)


class Document:
    def __init__(self, filename):
        """A container for a DOM-style SVG document.

        The `Document` class provides a simple interface to modify and analyze 
        the path elements in a DOM-style document.  The DOM-style document is 
        parsed into an ElementTree object (stored in the `tree` attribute).

        This class provides functions for extracting SVG data into Path objects.
        The Path output objects will be transformed based on their parent groups.
        
        Args:
            filename (str): The filename of the DOM-style object.
        """

        # remember location of original svg file
        if filename is not None and os.path.dirname(filename) == '':
            self.original_filename = os.path.join(os.getcwd(), filename)
        else:
            self.original_filename = filename

        if filename is not None:
            # parse svg to ElementTree object
            self.tree = etree.parse(filename)
        else:
            self.tree = etree.ElementTree(Element('svg'))

        self.root = self.tree.getroot()

    def flatten_all_paths(self, group_filter=lambda x: True,
                          path_filter=lambda x: True,
                          path_conversions=CONVERSIONS):
        """Forward the tree of this document into the more general
        flatten_all_paths function and return the result."""
        return flatten_all_paths(self.tree.getroot(), group_filter,
                                 path_filter, path_conversions)

    def flatten_group(self, group, recursive=True, group_filter=lambda x: True,
                      path_filter=lambda x: True, path_conversions=CONVERSIONS):
        if all(isinstance(s, str) for s in group):
            # If we're given a list of strings, assume it represents a
            # nested sequence
            group = self.get_or_add_group(group)
        elif not isinstance(group, Element):
            raise TypeError(
                'Must provide a list of strings that represent a nested '
                'group name, or provide an xml.etree.Element object. '
                'Instead you provided {0}'.format(group))

        return flatten_group(group, self.tree.getroot(), recursive,
                             group_filter, path_filter, path_conversions)

    def add_path(self, path, attribs=None, group=None):
        """Add a new path to the SVG."""

        # If not given a parent, assume that the path does not have a group
        if group is None:
            group = self.tree.getroot()

        # If given a list of strings (one or more), assume it represents
        # a sequence of nested group names
        elif all(isinstance(elem, str) for elem in group):
            group = self.get_or_add_group(group)

        elif not isinstance(group, Element):
            raise TypeError(
                'Must provide a list of strings or an xml.etree.Element '
                'object. Instead you provided {0}'.format(group))

        else:
            # Make sure that the group belongs to this Document object
            if not self.contains_group(group):
                warnings.warn('The requested group does not belong to '
                              'this Document')

        # TODO: It might be better to use duck-typing here with a try-except
        if isinstance(path, Path):
            path_svg = path.d()
        elif is_path_segment(path):
            path_svg = Path(path).d()
        elif isinstance(path, str):
            # Assume this is a valid d-string.
            # TODO: Should we sanity check the input string?
            path_svg = path
        else:
            raise TypeError(
                'Must provide a Path, a path segment type, or a valid '
                'SVG path d-string. Instead you provided {0}'.format(path))

        if attribs is None:
            attribs = {}
        else:
            attribs = attribs.copy()

        attribs['d'] = path_svg

        return SubElement(group, 'path', attribs)

    def contains_group(self, group):
        return any(group is owned for owned in self.tree.iter())

    def get_or_add_group(self, nested_names, name_attr='id'):
        """Get a group from the tree, or add a new one with the given
        name structure.

        `nested_names` is a list of strings which represent group names.
        Each group name will be nested inside of the previous group name.

        `name_attr` is the group attribute that is being used to
        represent the group's name. Default is 'id', but some SVGs may
        contain custom name labels, like 'inkscape:label'.

        Returns the requested group. If the requested group did not
        exist, this function will create it, as well as all parent
        groups that it requires. All created groups will be left with
        blank attributes.

        """
        group = self.tree.getroot()
        # Drill down through the names until we find the desired group
        while len(nested_names):
            prev_group = group
            next_name = nested_names.pop(0)
            for elem in group.iterfind(SVG_GROUP_TAG, SVG_NAMESPACE):
                if elem.get(name_attr) == next_name:
                    group = elem
                    break

            if prev_group is group:
                # The group we're looking for does not exist, so let's
                # create the group structure
                nested_names.insert(0, next_name)

                while nested_names:
                    next_name = nested_names.pop(0)
                    group = self.add_group({'id': next_name}, group)
                # Now nested_names will be empty, so the topmost
                # while-loop will end
        return group

    def add_group(self, group_attribs=None, parent=None):
        """Add an empty group element to the SVG."""
        if parent is None:
            parent = self.tree.getroot()
        elif not self.contains_group(parent):
            warnings.warn('The requested group {0} does not belong to '
                          'this Document'.format(parent))

        if group_attribs is None:
            group_attribs = {}
        else:
            group_attribs = group_attribs.copy()

        return SubElement(parent, '{{{0}}}g'.format(
            SVG_NAMESPACE['svg']), group_attribs)

    def save(self, filename):
        with open(filename, 'w') as output_svg:
            output_svg.write(etree.tostring(self.tree.getroot()))

    def display(self, filename=None):
        """Displays/opens the doc using the OS's default application."""

        if filename is None:
            filename = self.original_filename

        # write to a (by default temporary) file
        with open(filename, 'w') as output_svg:
            output_svg.write(etree.tostring(self.tree.getroot()))

        open_in_browser(filename)
