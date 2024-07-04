"""Tests related to SVG groups.

To run these tests, you can use (from root svgpathtools directory):
$ python -m unittest test.test_groups.TestGroups.test_group_flatten
"""
from __future__ import division, absolute_import, print_function
import unittest
from svgpathtools import Document, SVG_NAMESPACE, parse_path, Line, Arc
from os.path import join, dirname
import numpy as np


# When an assert fails, show the full error message, don't truncate it.
unittest.util._MAX_LENGTH = 999999999


def get_desired_path(name, paths):
    return next(p for p in paths
                if p.element.get('{some://testuri}name') == name)


class TestGroups(unittest.TestCase):

    def check_values(self, v, z):
        # Check that the components of 2D vector v match the components
        # of complex number z
        self.assertAlmostEqual(v[0], z.real)
        self.assertAlmostEqual(v[1], z.imag)

    def check_line(self, tf, v_s_vals, v_e_relative_vals, name, paths):
        # Check that the endpoints of the line have been correctly transformed.
        # * tf is the transform that should have been applied.
        # * v_s_vals is a 2D list of the values of the line's start point
        # * v_e_relative_vals is a 2D list of the values of the line's
        #    end point relative to the start point
        # * name is the path name (value of the test:name attribute in
        #    the SVG document)
        # * paths is the output of doc.paths()
        v_s_vals.append(1.0)
        v_e_relative_vals.append(0.0)
        v_s = np.array(v_s_vals)
        v_e = v_s + v_e_relative_vals

        actual = get_desired_path(name, paths)

        self.check_values(tf.dot(v_s), actual.start)
        self.check_values(tf.dot(v_e), actual.end)

    def test_group_transform(self):
        # The input svg has a group transform of "scale(1,-1)", which
        # can mess with Arc sweeps.
        doc = Document(join(dirname(__file__), 'negative-scale.svg'))
        path = doc.paths()[0]
        self.assertEqual(path[0], Line(start=-10j, end=-80j))
        self.assertEqual(path[1], Arc(start=-80j, radius=(30+30j), rotation=0.0, large_arc=True, sweep=True, end=-140j))
        self.assertEqual(path[2], Arc(start=-140j, radius=(20+20j), rotation=0.0, large_arc=False, sweep=False, end=-100j))
        self.assertEqual(path[3], Line(start=-100j, end=(100-100j)))
        self.assertEqual(path[4], Arc(start=(100-100j), radius=(20+20j), rotation=0.0, large_arc=True, sweep=False, end=(100-140j)))
        self.assertEqual(path[5], Arc(start=(100-140j), radius=(30+30j), rotation=0.0, large_arc=False, sweep=True, end=(100-80j)))
        self.assertEqual(path[6], Line(start=(100-80j), end=(100-10j)))
        self.assertEqual(path[7], Arc(start=(100-10j), radius=(10+10j), rotation=0.0, large_arc=False, sweep=True, end=(90+0j)))
        self.assertEqual(path[8], Line(start=(90+0j), end=(10+0j)))
        self.assertEqual(path[9], Arc(start=(10+0j), radius=(10+10j), rotation=0.0, large_arc=False, sweep=True, end=-10j))

    def test_group_flatten(self):
        # Test the Document.paths() function against the
        # groups.svg test file.
        # There are 12 paths in that file, with various levels of being
        # nested inside of group transforms.
        # The check_line function is used to reduce the boilerplate,
        # since all the tests are very similar.
        # This test covers each of the different types of transforms
        # that are specified by the SVG standard.
        doc = Document(join(dirname(__file__), 'groups.svg'))

        result = doc.paths()
        self.assertEqual(12, len(result))

        tf_matrix_group = np.array([[1.5, 0.0, -40.0],
                                    [0.0, 0.5, 20.0],
                                    [0.0, 0.0, 1.0]])

        self.check_line(tf_matrix_group,
                        [183, 183], [0.0, -50],
                        'path00', result)

        tf_scale_group = np.array([[1.25, 0.0, 0.0],
                                   [0.0, 1.25, 0.0],
                                   [0.0, 0.0, 1.0]])

        self.check_line(tf_matrix_group.dot(tf_scale_group),
                        [122, 320], [-50.0, 0.0],
                        'path01', result)

        self.check_line(tf_matrix_group.dot(tf_scale_group),
                        [150, 200], [-50, 25],
                        'path02', result)

        self.check_line(tf_matrix_group.dot(tf_scale_group),
                        [150, 200], [-50, 25],
                        'path03', result)

        tf_nested_translate_group = np.array([[1, 0, 20],
                                              [0, 1, 0],
                                              [0, 0, 1]])

        self.check_line(tf_matrix_group.dot(tf_scale_group
                                      ).dot(tf_nested_translate_group),
                        [150, 200], [-50, 25],
                        'path04', result)

        tf_nested_translate_xy_group = np.array([[1, 0, 20],
                                                 [0, 1, 30],
                                                 [0, 0, 1]])

        self.check_line(tf_matrix_group.dot(tf_scale_group
                                      ).dot(tf_nested_translate_xy_group),
                        [150, 200], [-50, 25],
                        'path05', result)

        tf_scale_xy_group = np.array([[0.5, 0, 0],
                                      [0, 1.5, 0.0],
                                      [0, 0, 1]])

        self.check_line(tf_matrix_group.dot(tf_scale_xy_group),
                        [122, 320], [-50, 0],
                        'path06', result)

        a_07 = 20.0*np.pi/180.0
        tf_rotate_group = np.array([[np.cos(a_07), -np.sin(a_07), 0],
                                     [np.sin(a_07),  np.cos(a_07), 0],
                                     [0, 0, 1]])

        self.check_line(tf_matrix_group.dot(tf_rotate_group),
                        [183, 183], [0, 30],
                        'path07', result)

        a_08 = 45.0*np.pi/180.0
        tf_rotate_xy_group_R = np.array([[np.cos(a_08), -np.sin(a_08), 0],
                                         [np.sin(a_08),  np.cos(a_08), 0],
                                         [0, 0, 1]])
        tf_rotate_xy_group_T = np.array([[1, 0, 183],
                                         [0, 1, 183],
                                         [0, 0, 1]])
        tf_rotate_xy_group = tf_rotate_xy_group_T.dot(
            tf_rotate_xy_group_R).dot(
            np.linalg.inv(tf_rotate_xy_group_T))

        self.check_line(tf_matrix_group.dot(tf_rotate_xy_group),
                        [183, 183], [0, 30],
                        'path08', result)

        a_09 = 5.0*np.pi/180.0
        tf_skew_x_group = np.array([[1, np.tan(a_09), 0],
                                    [0, 1, 0],
                                    [0, 0, 1]])

        self.check_line(tf_matrix_group.dot(tf_skew_x_group),
                        [183, 183], [40, 40],
                        'path09', result)

        a_10 = 5.0*np.pi/180.0
        tf_skew_y_group = np.array([[1, 0, 0],
                                    [np.tan(a_10), 1, 0],
                                    [0, 0, 1]])

        self.check_line(tf_matrix_group.dot(tf_skew_y_group),
                        [183, 183], [40, 40],
                        'path10', result)

        # This last test is for handling transforms that are defined as
        # attributes of a <path> element.
        a_11 = -40*np.pi/180.0
        tf_path11_R = np.array([[np.cos(a_11), -np.sin(a_11), 0],
                                 [np.sin(a_11),  np.cos(a_11), 0],
                                 [0, 0, 1]])
        tf_path11_T = np.array([[1, 0, 100],
                                [0, 1, 100],
                                [0, 0, 1]])
        tf_path11 = tf_path11_T.dot(tf_path11_R).dot(np.linalg.inv(tf_path11_T))

        self.check_line(tf_matrix_group.dot(tf_skew_y_group).dot(tf_path11),
                        [180, 20], [-70, 80],
                        'path11', result)

    def check_group_count(self, doc, expected_count):
        count = 0
        for _ in doc.tree.getroot().iter('{{{0}}}g'.format(SVG_NAMESPACE['svg'])):
            count += 1

        self.assertEqual(expected_count, count)

    def test_nested_group(self):
        # A bug in the flattened_paths_from_group() implementation made it so that only top-level
        # groups could have their paths flattened. This is a regression test to make
        # sure that when a nested group is requested, its paths can also be flattened.
        doc = Document(join(dirname(__file__), 'groups.svg'))
        result = doc.paths_from_group(['matrix group', 'scale group'])
        self.assertEqual(len(result), 5)

    def test_add_group(self):
        # Test `Document.add_group()` function and related Document functions.
        doc = Document(None)
        self.check_group_count(doc, 0)

        base_group = doc.add_group()
        base_group.set('id', 'base_group')
        self.assertTrue(doc.contains_group(base_group))
        self.check_group_count(doc, 1)

        child_group = doc.add_group(parent=base_group)
        child_group.set('id', 'child_group')
        self.assertTrue(doc.contains_group(child_group))
        self.check_group_count(doc, 2)

        grandchild_group = doc.add_group(parent=child_group)
        grandchild_group.set('id', 'grandchild_group')
        self.assertTrue(doc.contains_group(grandchild_group))
        self.check_group_count(doc, 3)

        sibling_group = doc.add_group(parent=base_group)
        sibling_group.set('id', 'sibling_group')
        self.assertTrue(doc.contains_group(sibling_group))
        self.check_group_count(doc, 4)

        # Test that we can retrieve each new group from the document
        self.assertEqual(base_group, doc.get_or_add_group(['base_group']))
        self.assertEqual(child_group, doc.get_or_add_group(
            ['base_group', 'child_group']))
        self.assertEqual(grandchild_group, doc.get_or_add_group(
            ['base_group', 'child_group', 'grandchild_group']))
        self.assertEqual(sibling_group, doc.get_or_add_group(
            ['base_group', 'sibling_group']))

        # Create a new nested group
        new_child = doc.get_or_add_group(
            ['base_group', 'new_parent', 'new_child'])
        self.check_group_count(doc, 6)
        self.assertEqual(new_child, doc.get_or_add_group(
            ['base_group', 'new_parent', 'new_child']))

        new_leaf = doc.get_or_add_group(
            ['base_group', 'new_parent', 'new_child', 'new_leaf'])
        self.assertEqual(new_leaf, doc.get_or_add_group([
            'base_group', 'new_parent', 'new_child', 'new_leaf']))
        self.check_group_count(doc, 7)

        path_d = ('M 206.071,858.413 L 206.071,-2.02031 '
                  'C -50.738,-81.1481 -20.364,-105.871 52.5279,-101.015 '
                  'L 103.036,0 '
                  'L 0,111.117')

        svg_path = doc.add_path(path_d, group=new_leaf)
        self.assertEqual(path_d, svg_path.get('d'))

        path = parse_path(path_d)
        svg_path = doc.add_path(path, group=new_leaf)
        self.assertEqual(path_d, svg_path.get('d'))

        # Test that paths are added to the correct group
        new_sibling = doc.get_or_add_group(
            ['base_group', 'new_parent', 'new_sibling'])
        doc.add_path(path, group=new_sibling)
        self.assertEqual(len(new_sibling), 1)
        self.assertEqual(path_d, new_sibling[0].get('d'))
