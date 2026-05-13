import unittest
import sys
import os

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _root not in sys.path:
    sys.path.insert(0, _root)

from tools.SegmentTree import SegmentTree


class TestSegmentTreeInit(unittest.TestCase):

    def test_create_with_coordinates(self):
        st = SegmentTree(0, 10)
        self.assertEqual(st.x1, 0)
        self.assertEqual(st.x2, 10)

    def test_create_with_none(self):
        st = SegmentTree()
        self.assertIsNone(st.x1)
        self.assertIsNone(st.x2)

    def test_create_swaps_if_x1_gt_x2(self):
        st = SegmentTree(10, 0)
        self.assertEqual(st.x1, 0)
        self.assertEqual(st.x2, 10)

    def test_create_with_name(self):
        st = SegmentTree(0, 10, name="root")
        self.assertEqual(st.name, "root")

    def test_children_initially_empty(self):
        st = SegmentTree(0, 10)
        self.assertEqual(st.children, [])


class TestSegmentTreeInsert(unittest.TestCase):

    def test_insert_child(self):
        root = SegmentTree(0, 100)
        child = root.insert(10, 50, name="child1")
        self.assertEqual(child.x1, 10)
        self.assertEqual(child.x2, 50)
        self.assertEqual(len(root.children), 1)

    def test_insert_multiple_children(self):
        root = SegmentTree(0, 100)
        root.insert(10, 30)
        root.insert(40, 60)
        root.insert(70, 90)
        self.assertEqual(len(root.children), 3)

    def test_insert_nested_child(self):
        root = SegmentTree(0, 100)
        root.insert(10, 50)
        root.insert(20, 30)
        # 20-30 should be nested inside 10-50
        self.assertEqual(len(root.children), 1)
        self.assertEqual(len(root.children[0].children), 1)
        self.assertEqual(root.children[0].children[0].x1, 20)

    def test_insert_duplicate_merges(self):
        root = SegmentTree(0, 100)
        root.insert(10, 50, name="first")
        root.insert(10, 50, name="second")
        self.assertEqual(len(root.children), 1)
        self.assertIn("first", root.children[0].name)
        self.assertIn("second", root.children[0].name)

    def test_insert_swaps_x1_x2(self):
        root = SegmentTree(0, 100)
        child = root.insert(50, 10, name="swapped")
        self.assertEqual(child.x1, 10)
        self.assertEqual(child.x2, 50)

    def test_insert_with_refered_object(self):
        root = SegmentTree(0, 100)
        child = root.insert(10, 50, referedObject="gene1")
        self.assertEqual(child.referedObject, ["gene1"])


class TestSegmentTreeIntersect(unittest.TestCase):

    def test_intersect_point(self):
        root = SegmentTree(0, 100)
        root.insert(10, 50, name="seg1")
        root.insert(60, 90, name="seg2")
        result = root.intersect(25)
        names = [r.name for r in result]
        self.assertTrue(any("seg1" in n for n in names))

    def test_intersect_range(self):
        root = SegmentTree(0, 100)
        root.insert(10, 50, name="seg1")
        root.insert(60, 90, name="seg2")
        result = root.intersect(45, 65)
        names = [r.name for r in result]
        self.assertTrue(any("seg1" in n for n in names))
        self.assertTrue(any("seg2" in n for n in names))

    def test_intersect_no_match(self):
        root = SegmentTree(0, 100)
        root.insert(10, 20)
        root.insert(30, 40)
        # Point between segments
        result = root.intersect(25)
        self.assertEqual(result, [])

    def test_intersect_empty_tree(self):
        root = SegmentTree(0, 100)
        result = root.intersect(50)
        self.assertEqual(result, [])


class TestSegmentTreeBounds(unittest.TestCase):

    def test_getX1_with_coordinates(self):
        st = SegmentTree(5, 10)
        self.assertEqual(st.getX1(), 5)

    def test_getX2_with_coordinates(self):
        st = SegmentTree(5, 10)
        self.assertEqual(st.getX2(), 10)

    def test_getX1_from_children(self):
        root = SegmentTree()
        root.insert(10, 50)
        root.insert(60, 90)
        self.assertEqual(root.getX1(), 10)

    def test_getX2_from_children(self):
        root = SegmentTree()
        root.insert(10, 50)
        root.insert(60, 90)
        self.assertEqual(root.getX2(), 90)


class TestSegmentTreeIndexedLength(unittest.TestCase):

    def test_length_with_coordinates(self):
        st = SegmentTree(0, 100)
        self.assertEqual(st.getIndexedLength(), 100)

    def test_length_no_children_no_coords(self):
        st = SegmentTree()
        self.assertEqual(st.getIndexedLength(), 0)

    def test_length_from_children(self):
        root = SegmentTree()
        root.insert(0, 10)
        root.insert(20, 30)
        self.assertEqual(root.getIndexedLength(), 20)


class TestSegmentTreeFirstLevel(unittest.TestCase):

    def test_first_level(self):
        root = SegmentTree(0, 100)
        root.insert(10, 30)
        root.insert(50, 70)
        result = root.getFirstLevel()
        self.assertEqual(result, [(10, 30), (50, 70)])

    def test_first_level_no_children(self):
        st = SegmentTree(0, 100)
        self.assertEqual(st.getFirstLevel(), [(0, 100)])

    def test_first_level_none_root(self):
        st = SegmentTree()
        self.assertIsNone(st.getFirstLevel())


class TestSegmentTreeEmptyChildren(unittest.TestCase):

    def test_empty_children(self):
        root = SegmentTree(0, 100)
        root.insert(10, 50)
        root.insert(60, 90)
        self.assertEqual(len(root.children), 2)
        root.emptyChildren()
        self.assertEqual(len(root.children), 0)


class TestSegmentTreeDichotomicSearch(unittest.TestCase):
    """Test that integer division fix in __dichotomicSearch works correctly."""

    def test_search_with_many_children(self):
        root = SegmentTree(0, 1000)
        for i in range(0, 100, 10):
            root.insert(i, i + 5, name=f"seg_{i}")
        # If dichotomicSearch used float division, intersect would fail with TypeError
        result = root.intersect(22)
        names = [r.name for r in result]
        self.assertTrue(any("seg_20" in n for n in names))


if __name__ == "__main__":
    unittest.main()
