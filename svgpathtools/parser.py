"""This submodule contains the path_parse() function used to convert SVG path
element d-strings into svgpathtools Path objects.
Note: This file was taken (nearly) as is from the svg.path module (v 2.0)."""

# External dependencies
from __future__ import division, absolute_import, print_function
import re
import numpy as np
import warnings

# Internal dependencies
from .path import Path, Line, QuadraticBezier, CubicBezier, Arc

# To maintain forward/backward compatibility
try:
    str = basestring
except NameError:
    pass

COMMANDS = set('MmZzLlHhVvCcSsQqTtAa')
UPPERCASE = set('MZLHVCSQTA')

COMMAND_RE = re.compile("([MmZzLlHhVvCcSsQqTtAa])")
FLOAT_RE = re.compile("[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?")


def _tokenize_path(pathdef):
    for x in COMMAND_RE.split(pathdef):
        if x in COMMANDS:
            yield x
        for token in FLOAT_RE.findall(x):
            yield token


def parse_path(pathdef, current_pos=0j, tree_element=None):
    # In the SVG specs, initial movetos are absolute, even if
    # specified as 'm'. This is the default behavior here as well.
    # But if you pass in a current_pos variable, the initial moveto
    # will be relative to that current_pos. This is useful.
    elements = list(_tokenize_path(pathdef))
    # Reverse for easy use of .pop()
    elements.reverse()

    if tree_element is None:
        segments = Path()
    else:
        segments = Path(tree_element=tree_element)

    start_pos = None
    command = None

    while elements:

        if elements[-1] in COMMANDS:
            # New command.
            last_command = command  # Used by S and T
            command = elements.pop()
            absolute = command in UPPERCASE
            command = command.upper()
        else:
            # If this element starts with numbers, it is an implicit command
            # and we don't change the command. Check that it's allowed:
            if command is None:
                raise ValueError("Unallowed implicit command in %s, position %s" % (
                    pathdef, len(pathdef.split()) - len(elements)))

        if command == 'M':
            # Moveto command.
            x = elements.pop()
            y = elements.pop()
            pos = float(x) + float(y) * 1j
            if absolute:
                current_pos = pos
            else:
                current_pos += pos

            # when M is called, reset start_pos
            # This behavior of Z is defined in svg spec:
            # http://www.w3.org/TR/SVG/paths.html#PathDataClosePathCommand
            start_pos = current_pos

            # Implicit moveto commands are treated as lineto commands.
            # So we set command to lineto here, in case there are
            # further implicit commands after this moveto.
            command = 'L'

        elif command == 'Z':
            # Close path
            if not (current_pos == start_pos):
                segments.append(Line(current_pos, start_pos))
            segments.closed = True
            current_pos = start_pos
            command = None

        elif command == 'L':
            x = elements.pop()
            y = elements.pop()
            pos = float(x) + float(y) * 1j
            if not absolute:
                pos += current_pos
            segments.append(Line(current_pos, pos))
            current_pos = pos

        elif command == 'H':
            x = elements.pop()
            pos = float(x) + current_pos.imag * 1j
            if not absolute:
                pos += current_pos.real
            segments.append(Line(current_pos, pos))
            current_pos = pos

        elif command == 'V':
            y = elements.pop()
            pos = current_pos.real + float(y) * 1j
            if not absolute:
                pos += current_pos.imag * 1j
            segments.append(Line(current_pos, pos))
            current_pos = pos

        elif command == 'C':
            control1 = float(elements.pop()) + float(elements.pop()) * 1j
            control2 = float(elements.pop()) + float(elements.pop()) * 1j
            end = float(elements.pop()) + float(elements.pop()) * 1j

            if not absolute:
                control1 += current_pos
                control2 += current_pos
                end += current_pos

            segments.append(CubicBezier(current_pos, control1, control2, end))
            current_pos = end

        elif command == 'S':
            # Smooth curve. First control point is the "reflection" of
            # the second control point in the previous path.

            if last_command not in 'CS':
                # If there is no previous command or if the previous command
                # was not an C, c, S or s, assume the first control point is
                # coincident with the current point.
                control1 = current_pos
            else:
                # The first control point is assumed to be the reflection of
                # the second control point on the previous command relative
                # to the current point.
                control1 = current_pos + current_pos - segments[-1].control2

            control2 = float(elements.pop()) + float(elements.pop()) * 1j
            end = float(elements.pop()) + float(elements.pop()) * 1j

            if not absolute:
                control2 += current_pos
                end += current_pos

            segments.append(CubicBezier(current_pos, control1, control2, end))
            current_pos = end

        elif command == 'Q':
            control = float(elements.pop()) + float(elements.pop()) * 1j
            end = float(elements.pop()) + float(elements.pop()) * 1j

            if not absolute:
                control += current_pos
                end += current_pos

            segments.append(QuadraticBezier(current_pos, control, end))
            current_pos = end

        elif command == 'T':
            # Smooth curve. Control point is the "reflection" of
            # the second control point in the previous path.

            if last_command not in 'QT':
                # If there is no previous command or if the previous command
                # was not an Q, q, T or t, assume the first control point is
                # coincident with the current point.
                control = current_pos
            else:
                # The control point is assumed to be the reflection of
                # the control point on the previous command relative
                # to the current point.
                control = current_pos + current_pos - segments[-1].control

            end = float(elements.pop()) + float(elements.pop()) * 1j

            if not absolute:
                end += current_pos

            segments.append(QuadraticBezier(current_pos, control, end))
            current_pos = end

        elif command == 'A':
            radius = float(elements.pop()) + float(elements.pop()) * 1j
            rotation = float(elements.pop())
            arc = float(elements.pop())
            sweep = float(elements.pop())
            end = float(elements.pop()) + float(elements.pop()) * 1j

            if not absolute:
                end += current_pos

            segments.append(Arc(current_pos, radius, rotation, arc, sweep, end))
            current_pos = end

    return segments


def _check_num_parsed_values(values, allowed):
    if not any(num == len(values) for num in allowed):
        if len(allowed) > 1:
            warnings.warn('Expected one of the following number of values {0}, but found {1} values instead: {2}'
                          .format(allowed, len(values), values))
        elif allowed[0] != 1:
            warnings.warn('Expected {0} values, found {1}: {2}'.format(allowed[0], len(values), values))
        else:
            warnings.warn('Expected 1 value, found {0}: {1}'.format(len(values), values))
        return False
    return True


def _parse_transform_substr(transform_substr):

    type_str, value_str = transform_substr.split('(')
    value_str = value_str.replace(',', ' ')
    values = list(map(float, filter(None, value_str.split(' '))))

    transform = np.identity(3)
    if 'matrix' in type_str:
        if not _check_num_parsed_values(values, [6]):
            return transform

        transform[0:2, 0:3] = np.array([values[0:6:2], values[1:6:2]])

    elif 'translate' in transform_substr:
        if not _check_num_parsed_values(values, [1, 2]):
            return transform

        transform[0, 2] = values[0]
        if len(values) > 1:
            transform[1, 2] = values[1]

    elif 'scale' in transform_substr:
        if not _check_num_parsed_values(values, [1, 2]):
            return transform

        x_scale = values[0]
        y_scale = values[1] if (len(values) > 1) else x_scale
        transform[0, 0] = x_scale
        transform[1, 1] = y_scale

    elif 'rotate' in transform_substr:
        if not _check_num_parsed_values(values, [1, 3]):
            return transform

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
        if not _check_num_parsed_values(values, [1]):
            return transform

        transform[0, 1] = np.tan(values[0] * np.pi / 180.0)

    elif 'skewY' in transform_substr:
        if not _check_num_parsed_values(values, [1]):
            return transform

        transform[1, 0] = np.tan(values[0] * np.pi / 180.0)
    else:
        # Return an identity matrix if the type of transform is unknown, and warn the user
        warnings.warn('Unknown SVG transform type: {0}'.format(type_str))

    return transform


def parse_transform(transform_str):
    """Converts a valid SVG transformation string into a 3x3 matrix.
    If the string is empty or null, this returns a 3x3 identity matrix"""
    if not transform_str:
        return np.identity(3)
    elif not isinstance(transform_str, str):
        raise TypeError('Must provide a string to parse')

    total_transform = np.identity(3)
    transform_substrs = transform_str.split(')')[:-1]  # Skip the last element, because it should be empty
    for substr in transform_substrs:
        total_transform = total_transform.dot(_parse_transform_substr(substr))

    return total_transform
