import unittest
import os
import sys
import tempfile

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _root not in sys.path:
    sys.path.insert(0, _root)

from tools.parsers.FastaTools import FastaFile


class TestFastaFileCreation(unittest.TestCase):

    def test_create_empty(self):
        f = FastaFile()
        self.assertEqual(len(f), 0)

    def test_add_entry(self):
        f = FastaFile()
        f.add(">seq1", "ATCGATCG")
        self.assertEqual(len(f), 1)

    def test_add_entry_without_gt(self):
        f = FastaFile()
        f.add("seq1", "ATCGATCG")
        entry = f.get(0)
        self.assertTrue(entry[0].startswith(">"))

    def test_add_multiple(self):
        f = FastaFile()
        f.add(">seq1", "AAAA")
        f.add(">seq2", "CCCC")
        f.add(">seq3", "GGGG")
        self.assertEqual(len(f), 3)


class TestFastaFileGet(unittest.TestCase):

    def test_get_entry(self):
        f = FastaFile()
        f.add(">seq1", "ATCG")
        entry = f.get(0)
        self.assertEqual(entry[0], ">seq1")
        self.assertEqual(entry[1], "ATCG")

    def test_getitem(self):
        f = FastaFile()
        f.add(">seq1", "ATCG")
        entry = f[0]
        self.assertEqual(entry[1], "ATCG")


class TestFastaFileParseStr(unittest.TestCase):

    def test_parse_simple(self):
        f = FastaFile()
        f.parseStr(">seq1\nATCG\n>seq2\nGGCC\n")
        self.assertEqual(len(f), 2)
        entry0 = f.get(0)
        self.assertEqual(entry0[0], ">seq1")
        self.assertEqual(entry0[1], "ATCG")

    def test_parse_multiline_sequence(self):
        f = FastaFile()
        f.parseStr(">seq1\nATCG\nGGCC\n")
        entry = f.get(0)
        self.assertEqual(entry[1], "ATCGGGCC")


class TestFastaFileSaveAndParse(unittest.TestCase):

    def setUp(self):
        self.tmpfile = tempfile.mktemp(suffix=".fasta")

    def tearDown(self):
        if os.path.exists(self.tmpfile):
            os.remove(self.tmpfile)

    def test_save_and_parse_roundtrip(self):
        f = FastaFile()
        f.add(">gene1", "ATGATGATG")
        f.add(">gene2", "CCCGGGAAA")
        f.save(self.tmpfile)

        f2 = FastaFile(self.tmpfile)
        self.assertEqual(len(f2), 2)
        self.assertEqual(f2.get(0)[0], ">gene1")
        self.assertEqual(f2.get(0)[1], "ATGATGATG")
        self.assertEqual(f2.get(1)[0], ">gene2")
        self.assertEqual(f2.get(1)[1], "CCCGGGAAA")


class TestFastaFileIteration(unittest.TestCase):

    def test_iteration(self):
        f = FastaFile()
        f.add(">s1", "AAA")
        f.add(">s2", "CCC")
        f.add(">s3", "GGG")
        entries = list(f)
        self.assertEqual(len(entries), 3)
        self.assertEqual(entries[0][1], "AAA")
        self.assertEqual(entries[2][1], "GGG")


class TestFastaFileToStr(unittest.TestCase):

    def test_to_str(self):
        f = FastaFile()
        f.add(">seq1", "ATCG")
        s = f.toStr()
        self.assertIn(">seq1", s)
        self.assertIn("ATCG", s)


class TestFastaFileSetItem(unittest.TestCase):

    def test_setitem(self):
        f = FastaFile()
        f.add(">seq1", "ATCG")
        f[0] = (">newseq", "GGGG")
        self.assertEqual(f[0][0], ">newseq")
        self.assertEqual(f[0][1], "GGGG")

    def test_setitem_wrong_length_raises(self):
        f = FastaFile()
        f.add(">seq1", "ATCG")
        with self.assertRaises(TypeError):
            f[0] = ("only_one_element",)


class TestFastaFileReset(unittest.TestCase):

    def test_reset(self):
        f = FastaFile()
        f.add(">s1", "AAA")
        f.add(">s2", "CCC")
        self.assertEqual(len(f), 2)
        f.reset()
        self.assertEqual(len(f), 0)


if __name__ == "__main__":
    unittest.main()
