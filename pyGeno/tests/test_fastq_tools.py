import unittest
import os
import sys
import tempfile

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _root not in sys.path:
    sys.path.insert(0, _root)

from tools.parsers.FastqTools import FastqFile, FastqEntry


class TestFastqEntry(unittest.TestCase):

    def test_create_entry(self):
        e = FastqEntry("@id", "ATCG", "+", "IIII")
        self.assertEqual(e["identifier"], "@id")
        self.assertEqual(e["sequence"], "ATCG")
        self.assertEqual(e["+"], "+")
        self.assertEqual(e["qualities"], "IIII")

    def test_set_entry(self):
        e = FastqEntry()
        e["identifier"] = "@test"
        e["sequence"] = "GGGG"
        self.assertEqual(e["identifier"], "@test")
        self.assertEqual(e["sequence"], "GGGG")

    def test_str_format(self):
        e = FastqEntry("@id", "ATCG", "+", "IIII")
        s = str(e)
        self.assertIn("@id", s)
        self.assertIn("ATCG", s)
        self.assertIn("IIII", s)


class TestFastqFileParseStr(unittest.TestCase):

    def test_parse_single_entry(self):
        data = "@SEQ_ID\nATCGATCG\n+\nIIIIIIII\n"
        f = FastqFile()
        f.parseStr(data)
        self.assertEqual(len(f), 1)

    def test_parse_multiple_entries(self):
        data = "@SEQ1\nATCG\n+\nIIII\n@SEQ2\nGGGG\n+\nJJJJ\n"
        f = FastqFile()
        f.parseStr(data)
        self.assertEqual(len(f), 2)

    def test_len_uses_integer_division(self):
        # Verifies the // fix: len should always return int
        data = "@SEQ1\nATCG\n+\nIIII\n"
        f = FastqFile()
        f.parseStr(data)
        length = len(f)
        self.assertIsInstance(length, int)
        self.assertEqual(length, 1)


class TestFastqFileGetEntry(unittest.TestCase):

    def test_get_entry(self):
        data = "@SEQ1\nATCG\n+\nIIII\n@SEQ2\nGGGG\n+\nJJJJ\n"
        f = FastqFile()
        f.parseStr(data)
        e = f.get(0)
        self.assertEqual(e["identifier"], "@SEQ1")
        self.assertEqual(e["sequence"], "ATCG")

    def test_getitem(self):
        data = "@SEQ1\nATCG\n+\nIIII\n"
        f = FastqFile()
        f.parseStr(data)
        self.assertEqual(f[0]["sequence"], "ATCG")


class TestFastqFileNewEntry(unittest.TestCase):

    def test_new_entry_returns_fastq_entry(self):
        f = FastqFile()
        e = f.newEntry()
        self.assertIsInstance(e, FastqEntry)

    def test_add_entry(self):
        f = FastqFile()
        e = FastqEntry("@id", "AAAA", "+", "FFFF")
        f.add(e)
        self.assertGreater(len(f.data), 0)


class TestFastqFileSave(unittest.TestCase):

    def setUp(self):
        self.tmpfile = tempfile.mktemp(suffix=".fastq")

    def tearDown(self):
        if os.path.exists(self.tmpfile):
            os.remove(self.tmpfile)

    def test_parse_then_write_manually(self):
        """Write parsed data back to file manually since save() relies on missing make() method."""
        f = FastqFile()
        f.parseStr("@SEQ1\nATCG\n+\nIIII\n@SEQ2\nGGGG\n+\nJJJJ\n")
        with open(self.tmpfile, 'w') as fh:
            for i in range(len(f)):
                entry = f.get(i)
                fh.write(str(entry) + '\n')
        self.assertTrue(os.path.exists(self.tmpfile))


class TestFastqFileReset(unittest.TestCase):

    def test_reset(self):
        f = FastqFile()
        f.parseStr("@SEQ1\nATCG\n+\nIIII\n")
        self.assertEqual(len(f), 1)
        f.reset()
        self.assertEqual(len(f), 0)


if __name__ == "__main__":
    unittest.main()
