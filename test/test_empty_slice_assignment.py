import unittest
from svgpathtools import Line, Path


class TestEmptySliceAssignment(unittest.TestCase):
    def test_clearing_path_by_slice_resets_endpoints_and_length(self):
        path = Path(Line(0j, 1j), Line(1j, 2j))
        self.assertEqual(path.length(), 2)
        path[:] = []
        self.assertEqual(len(path), 0)
        self.assertIsNone(path.start)
        self.assertIsNone(path.end)
        self.assertEqual(path.length(), 0)
        path.append(Line(3j, 5j))
        self.assertEqual(path.start, 3j)
        self.assertEqual(path.end, 5j)
        self.assertEqual(path.length(), 2)

    def test_empty_assignment_to_empty_path_is_valid(self):
        path = Path()
        path[:] = []
        self.assertEqual(path, Path())
        self.assertIsNone(path.start)
        self.assertIsNone(path.end)

    def test_partial_and_nonempty_slice_assignments_still_work(self):
        path = Path(Line(0j, 1j), Line(1j, 2j))
        path[:1] = []
        self.assertEqual(path.start, 1j)
        self.assertEqual(path.end, 2j)
        path[:] = [Line(4j, 7j)]
        self.assertEqual(path.length(), 3)
        self.assertEqual(path.start, 4j)
        self.assertEqual(path.end, 7j)
