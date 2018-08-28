"""(Experimental) replacement for import/export functionality SAX

"""

# External dependencies
from __future__ import division, absolute_import, print_function
import os
import xml.sax
from xml.etree.ElementTree import Element, ElementTree, SubElement

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

NAME_SVG = "svg"
ATTR_VERSION = "version"
VALUE_SVG_VERSION = "1.1"
ATTR_XMLNS = "xmlns"
VALUE_XMLNS = "http://www.w3.org/2000/svg"
ATTR_XMLNS_LINK = "xmlns:xlink"
VALUE_XLINK = "http://www.w3.org/1999/xlink"
ATTR_XMLNS_EV = "xmlns:ev"
VALUE_XMLNS_EV = "http://www.w3.org/2001/xml-events"
ATTR_WIDTH = "width"
ATTR_HEIGHT = "height"
ATTR_VIEWBOX = "viewBox"
NAME_PATH = "path"
ATTR_DATA = "d"
ATTR_FILL = "fill"
ATTR_STROKE = "stroke"
ATTR_STROKE_WIDTH = "stroke-width"
VALUE_NONE = "none"


class SvgIoSaxHandler(xml.sax.ContentHandler):
    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.root_values = {}
        self.tree = []
        self.stack = []
        self.values = {}
        self.matrix = np.identity(3)

    def push_stack(self):
        self.stack.append((self.values, self.matrix))
        self.matrix = np.identity(3).dot(self.matrix)  # copy of matrix
        current_values = self.values
        self.values = {}
        self.values.update(current_values)  # copy of dictionary

    def pop_stack(self):
        v = self.stack.pop()
        self.values = v[0]
        self.matrix = v[1]

    def startElement(self, name, attrs):
        self.push_stack()
        for (k, v) in attrs.items():
            if k is None:
                continue
            elif "style" == k:
                for equate in v.split(";"):
                    for p, q in equate:
                        self.values[p] = q
            elif "transform" == k:
                self.matrix.dot(parse_transform(v))
            else:
                self.values[k] = v
        if "svg" == name:
            current_values = self.values
            self.values = {}
            self.values.update(current_values)
            self.root_values = current_values
        elif "g" == name:
            return
        elif 'path' == name:
            self.tree.append((self.values, parse_path(self.values['d']), self.matrix))
        elif 'circle' == name:
            pass
        elif 'ellipse' == name:
            pass
        elif 'line' == name:
            pass
        elif 'polyline' == name:
            pass
        elif 'polygon' == name:
            pass
        elif 'rect' == name:
            pass

    def endElement(self, name):
        self.pop_stack()


class SaxDocument:
    def __init__(self, filename):
        """A container for a SAX SVG light tree objects document.

        This class provides functions for extracting SVG data into Path objects.

        Args:
            filename (str): The filename of the SVG file
        """

        # remember location of original svg file
        if filename is not None and os.path.dirname(filename) == '':
            self.original_filename = os.path.join(os.getcwd(), filename)
        else:
            self.original_filename = filename

        if filename is not None:
            parser = xml.sax.make_parser()
            handler = SvgIoSaxHandler()
            parser.setContentHandler(handler)
            parser.parse(open(filename, "r"))
            self.tree = handler.tree
            self.root_values = handler.root_values

    def generate_dom(self):
        root = Element(NAME_SVG)
        root.set(ATTR_VERSION, VALUE_SVG_VERSION)
        root.set(ATTR_XMLNS, VALUE_XMLNS)
        root.set(ATTR_XMLNS_LINK, VALUE_XLINK)
        root.set(ATTR_XMLNS_EV, VALUE_XMLNS_EV)
        width = self.root_values.get(ATTR_WIDTH, None)
        height = self.root_values.get(ATTR_HEIGHT, None)
        if width is not None:
            root.set(ATTR_WIDTH, width)
        if height is not None:
            root.set(ATTR_HEIGHT, height)
        viewbox = self.root_values.get(ATTR_VIEWBOX, None)
        if viewbox is not None:
            root.set(ATTR_VIEWBOX, viewbox)
        for element in self.tree:
            path = SubElement(root, NAME_PATH)
            path.set(ATTR_DATA, element[0][ATTR_DATA])
            path.set(ATTR_FILL, element[0][ATTR_FILL])
            path.set(ATTR_STROKE, element[0][ATTR_STROKE])
        return ElementTree(root)

    def save(self, filename):
        with open(filename, 'w') as output_svg:
            dom_tree = self.generate_dom()
            dom_tree.write(output_svg)

    def display(self, filename=None):
        """Displays/opens the doc using the OS's default application."""
        if filename is None:
            filename = 'display_temp.svg'
        self.save(filename)
        open_in_browser(filename)
