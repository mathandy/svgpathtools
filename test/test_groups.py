from __future__ import division, absolute_import, print_function
import unittest
from svgpathtools import *
from os.path import join, dirname
import numpy as np


def get_desired_path(name, paths):
    return next(p for p in paths if p.attributes['{some://testuri}name'] == name)


def column_vector(values):
    input = []
    for value in values:
        input.append([value])
    return np.matrix(input)


class TestGroups(unittest.TestCase):

    def check_values(self, v, z):
        self.assertAlmostEqual(v[0], z.real)
        self.assertAlmostEqual(v[1], z.imag)

    def check_path(self, tf, v_s_vals, v_e_relative_vals, name, paths):
        v_s_vals.append(1.0)
        v_e_relative_vals.append(0.0)
        v_s = column_vector(v_s_vals)
        v_e = v_s + column_vector(v_e_relative_vals)

        actual = get_desired_path(name, paths)

        self.check_values(tf.dot(v_s), actual.path.start)
        self.check_values(tf.dot(v_e), actual.path.end)

    def test_group_flatten(self):
        doc = Document(join(dirname(__file__), 'groups.svg'))

        result = doc.flatten_all_paths()
        self.assertEqual(11, len(result))

        tf_matrix_group = np.matrix([[1.5, 0.0, -40.0], [0.0, 0.5, 20.0], [0.0, 0.0, 1.0]])

        self.check_path(tf_matrix_group,
                        [183, 183], [0.0, -50],
                        'path00', result)

        tf_scale_group = np.matrix([[1.25, 0.0, 0.0], [0.0, 1.25, 0.0], [0.0, 0.0, 1.0]])

        self.check_path(tf_matrix_group.dot(tf_scale_group),
                        [122, 320], [-50.0, 0.0],
                        'path01', result)

        self.check_path(tf_matrix_group.dot(tf_scale_group),
                        [150, 200], [-50, 25],
                        'path02', result)

        self.check_path(tf_matrix_group.dot(tf_scale_group),
                        [150, 200], [-50, 25],
                        'path03', result)

        tf_nested_translate_group = np.matrix([[1, 0, 20], [0, 1, 0], [0, 0, 1]])

        self.check_path(tf_matrix_group.dot(tf_scale_group).dot(tf_nested_translate_group),
                        [150, 200], [-50, 25],
                        'path04', result)

        tf_nested_translate_xy_group = np.matrix([[1, 0, 20], [0, 1, 30], [0, 0, 1]])

        self.check_path(tf_matrix_group.dot(tf_scale_group).dot(tf_nested_translate_xy_group),
                        [150, 200], [-50, 25],
                        'path05', result)

        tf_scale_xy_group = np.matrix([[0.5, 0, 0], [0, 1.5, 0.0], [0, 0, 1]])

        self.check_path(tf_matrix_group.dot(tf_scale_xy_group),
                        [122, 320], [-50, 0],
                        'path06', result)

        a_07 = 20.0*np.pi/180.0
        tf_rotate_group = np.matrix([[np.cos(a_07), -np.sin(a_07), 0],
                                     [np.sin(a_07),  np.cos(a_07), 0],
                                     [0, 0, 1]])

        self.check_path(tf_matrix_group.dot(tf_rotate_group),
                        [183, 183], [0, 30],
                        'path07', result)

        a_08 = 45.0*np.pi/180.0
        tf_rotate_xy_group_R = np.matrix([[np.cos(a_08), -np.sin(a_08), 0],
                                          [np.sin(a_08),  np.cos(a_08), 0],
                                          [0, 0, 1]])
        tf_rotate_xy_group_T = np.matrix([[1, 0, 183], [0, 1, 183], [0, 0, 1]])
        tf_rotate_xy_group = tf_rotate_xy_group_T.dot(tf_rotate_xy_group_R).dot(np.linalg.inv(tf_rotate_xy_group_T))

        self.check_path(tf_matrix_group.dot(tf_rotate_xy_group),
                        [183, 183], [0, 30],
                        'path08', result)

        a_09 = 5.0*np.pi/180.0
        tf_skew_x_group = np.matrix([[1, np.tan(a_09), 0], [0, 1, 0], [0, 0, 1]])

        self.check_path(tf_matrix_group.dot(tf_skew_x_group),
                        [183, 183], [40, 40],
                        'path09', result)

        a_10 = 5.0*np.pi/180.0
        tf_skew_y_group = np.matrix([[1, 0, 0], [np.tan(a_10), 1, 0], [0, 0, 1]])

        self.check_path(tf_matrix_group.dot(tf_skew_y_group),
                        [183, 183], [40, 40],
                        'path10', result)

        # TODO: Add a test where a path element has a transform attribute
