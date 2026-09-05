"""This submodule contains the path_parse() function used to convert SVG path
element d-strings into svgpathtools Path objects.
Note: This file was taken (nearly) as is from the svg.path module (v 2.0)."""

# External dependencies
from __future__ import division, absolute_import, print_function
from typing import Optional, Sequence
import numpy as np
import warnings

# Internal dependencies
from .path import Path


def parse_path(pathdef, current_pos=0j, tree_element=None):
    return Path(pathdef, current_pos=current_pos, tree_element=tree_element)


def _check_num_parsed_values(values: Sequence[float], allowed: Sequence[int]) -> None:
    if not any(num == len(values) for num in allowed):
        if len(allowed) > 1:
            raise ValueError('Expected one of the following number of values {0}, but found {1} values instead: {2}'
                             .format(allowed, len(values), values))
        elif allowed[0] != 1:
            raise ValueError('Expected {0} values, found {1}: {2}'.format(allowed[0], len(values), values))
        else:
            raise ValueError('Expected 1 value, found {0}: {1}'.format(len(values), values))


def _parse_transform_substr(transform_substr: str) -> np.ndarray:
    """Converts a single SVG transform substring, `type(v1 v2 ...)`, into
    a 3x3 matrix.  Raises a ValueError on invalid transform syntax."""

    if transform_substr.count('(') != 1:
        raise ValueError('Invalid SVG transform substring: {0!r}'.format(transform_substr))

    type_str, value_str = transform_substr.split('(')
    try:
        values = [float(s) for s in value_str.replace(',', ' ').split()]
    except ValueError:
        raise ValueError('Invalid SVG transform substring: {0!r}'.format(transform_substr))

    transform = np.identity(3)
    if 'matrix' in type_str:
        _check_num_parsed_values(values, [6])

        transform[0:2, 0:3] = np.array([values[0:6:2], values[1:6:2]])

    elif 'translate' in transform_substr:
        _check_num_parsed_values(values, [1, 2])

        transform[0, 2] = values[0]
        if len(values) > 1:
            transform[1, 2] = values[1]

    elif 'scale' in transform_substr:
        _check_num_parsed_values(values, [1, 2])

        x_scale = values[0]
        y_scale = values[1] if (len(values) > 1) else x_scale
        transform[0, 0] = x_scale
        transform[1, 1] = y_scale

    elif 'rotate' in transform_substr:
        _check_num_parsed_values(values, [1, 3])

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

    elif 'skewX' in transform_substr:
        _check_num_parsed_values(values, [1])

        transform[0, 1] = np.tan(values[0] * np.pi / 180.0)

    elif 'skewY' in transform_substr:
        _check_num_parsed_values(values, [1])

        transform[1, 0] = np.tan(values[0] * np.pi / 180.0)
    else:
        raise ValueError('Unknown SVG transform type: {0}'.format(type_str))

    return transform


def parse_transform(transform_str: Optional[str], strict: bool = False) -> np.ndarray:
    """Converts a valid SVG transformation string into a 3x3 matrix.
    If the string is empty or null, this returns a 3x3 identity matrix.

    By default each invalid transform substring is skipped (i.e.
    contributes an identity matrix) with a warning.  If `strict` is
    true, a ValueError is raised on invalid transform syntax instead."""
    if not transform_str:
        return np.identity(3)
    elif not isinstance(transform_str, str):
        raise TypeError('Must provide a string to parse')

    transform_substrs = transform_str.split(')')
    # Anything after the last ')' (e.g. a stray 'matrix' with no
    # parentheses) is invalid syntax, not a transform to apply.
    trailing = transform_substrs.pop()
    if trailing.strip():
        if strict:
            raise ValueError('Invalid SVG transform substring: {0!r}'.format(trailing))
        warnings.warn('Skipping invalid SVG transform substring: {0!r}'.format(trailing))

    total_transform = np.identity(3)
    for substr in transform_substrs:
        try:
            total_transform = total_transform.dot(_parse_transform_substr(substr))
        except ValueError as e:
            if strict:
                raise
            warnings.warn('Skipping invalid SVG transform substring {0!r}: {1}'.format(substr, e))

    return total_transform
