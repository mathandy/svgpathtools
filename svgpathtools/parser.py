"""
Parse SVG path element d-strings into svgpathtools Path objects.

This submodule contains the parse_path() function.  Note: this file was
taken (nearly) as is from the svg.path module (v 2.0).
"""

# External dependencies
from __future__ import division, absolute_import, print_function
from typing import Optional, Sequence, TYPE_CHECKING
import math
import re
import numpy as np
import warnings

# Internal dependencies
from .path import Path

if TYPE_CHECKING:
    # Imported for type annotations only; no XML is parsed here.
    from xml.etree.ElementTree import Element

# The SVG number grammar. Stricter than float(), which also accepts
# e.g. '1_0', 'nan', and non-ASCII digits.  (The one deliberate
# looseness: '1.' is accepted, though the grammar technically wants a
# digit after the point.)
_NUMBER_RE = re.compile(r'[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?\Z')


class SVGSyntaxWarning(UserWarning):
    """Category for warnings about invalid SVG syntax handled leniently."""


def parse_path(pathdef: str, current_pos: complex = 0j,
               tree_element: Optional['Element'] = None) -> Path:
    """Convert an SVG path element d-string into a Path object."""
    return Path(pathdef, current_pos=current_pos, tree_element=tree_element)


def _check_num_parsed_values(values: Sequence[float], allowed: Sequence[int],
                             transform_substr: str) -> None:
    if not any(num == len(values) for num in allowed):
        if len(allowed) > 1:
            raise ValueError('Expected one of the following number of values {0}, '
                             'but found {1} values in {2!r}: {3}'
                             .format(allowed, len(values), transform_substr, values))
        elif allowed[0] != 1:
            raise ValueError('Expected {0} values in {1!r}, found {2}: {3}'
                             .format(allowed[0], transform_substr, len(values), values))
        else:
            raise ValueError('Expected 1 value in {0!r}, found {1}: {2}'
                             .format(transform_substr, len(values), values))


def _parse_transform_substr(transform_substr: str) -> np.ndarray:
    """
    Convert a single SVG transform substring into a 3x3 matrix.

    A well-formed transform substring has the form `type(v1 v2 ...)`.
    Raises a ValueError on invalid transform syntax.
    """
    if transform_substr.count('(') != 1:
        raise ValueError('Invalid SVG transform substring: {0!r}'.format(transform_substr))

    type_str, value_str = transform_substr.split('(')
    # Any leading commas/whitespace are the separator from the preceding
    # transform in the list, e.g. 'translate(1), rotate(30)'.
    type_str = type_str.strip(', \t\n\r')
    tokens = value_str.replace(',', ' ').split()
    if not all(_NUMBER_RE.match(t) for t in tokens):
        raise ValueError('Invalid SVG transform substring: {0!r}'.format(transform_substr))
    values = [float(t) for t in tokens]
    # A grammar-valid value can still overflow float, e.g. '1e999'.
    if not all(math.isfinite(v) for v in values):
        raise ValueError('Non-finite value in SVG transform substring: {0!r}'.format(transform_substr))

    transform = np.identity(3)
    if type_str == 'matrix':
        _check_num_parsed_values(values, [6], transform_substr)

        transform[0:2, 0:3] = np.array([values[0:6:2], values[1:6:2]])

    elif type_str == 'translate':
        _check_num_parsed_values(values, [1, 2], transform_substr)

        transform[0, 2] = values[0]
        if len(values) > 1:
            transform[1, 2] = values[1]

    elif type_str == 'scale':
        _check_num_parsed_values(values, [1, 2], transform_substr)

        x_scale = values[0]
        y_scale = values[1] if (len(values) > 1) else x_scale
        transform[0, 0] = x_scale
        transform[1, 1] = y_scale

    elif type_str == 'rotate':
        _check_num_parsed_values(values, [1, 3], transform_substr)

        angle = values[0] * np.pi / 180.0
        if len(values) == 3:
            offset = values[1:3]
        else:
            offset = (0, 0)
        tf_offset = np.identity(3)
        tf_offset[0:2, 2:3] = np.array([[offset[0]], [offset[1]]])
        tf_rotate = np.identity(3)
        tf_rotate[0:2, 0:2] = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        tf_offset_neg = np.identity(3)
        tf_offset_neg[0:2, 2:3] = np.array([[-offset[0]], [-offset[1]]])

        transform = tf_offset.dot(tf_rotate).dot(tf_offset_neg)

    elif type_str == 'skewX':
        _check_num_parsed_values(values, [1], transform_substr)

        transform[0, 1] = np.tan(values[0] * np.pi / 180.0)

    elif type_str == 'skewY':
        _check_num_parsed_values(values, [1], transform_substr)

        transform[1, 0] = np.tan(values[0] * np.pi / 180.0)
    else:
        raise ValueError('Unknown SVG transform type: {0}'.format(type_str))

    return transform


def parse_transform(transform_str: Optional[str], strict: bool = False) -> np.ndarray:
    """
    Convert a valid SVG transformation string into a 3x3 matrix.

    If the string is empty, null, or 'none', this returns a 3x3
    identity matrix.  By default each invalid transform substring is
    skipped (i.e. contributes an identity matrix) with a warning.  If
    `strict` is true, a ValueError is raised on invalid transform
    syntax instead.

    The warnings use the SVGSyntaxWarning category, so they can be
    silenced or escalated selectively, e.g.
    `warnings.simplefilter('error', SVGSyntaxWarning)`.
    """
    if transform_str is None or transform_str == '':
        return np.identity(3)
    elif not isinstance(transform_str, str):
        raise TypeError('Must provide a string to parse')

    # 'none' is valid (SVG 2 / CSS transform syntax) and means no transform.
    if transform_str.strip() == 'none':
        return np.identity(3)

    transform_substrs = transform_str.split(')')
    # Anything after the last ')' (e.g. a stray 'matrix' with no
    # parentheses) is invalid syntax, not a transform to apply -- but a
    # trailing list-separator comma is harmless.
    trailing = transform_substrs.pop()
    if trailing.strip(', \t\n\r'):
        if strict:
            raise ValueError('Invalid SVG transform substring: {0!r}'.format(trailing))
        warnings.warn('Skipping invalid SVG transform substring: {0!r}'.format(trailing),
                      SVGSyntaxWarning)

    total_transform = np.identity(3)
    for substr in transform_substrs:
        try:
            total_transform = total_transform.dot(_parse_transform_substr(substr))
        except ValueError as e:
            if strict:
                raise
            warnings.warn('Skipping invalid SVG transform substring {0!r}: {1}'.format(substr, e),
                          SVGSyntaxWarning)

    return total_transform
