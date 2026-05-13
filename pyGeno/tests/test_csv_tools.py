import unittest
import os
import sys
import tempfile

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _root not in sys.path:
    sys.path.insert(0, _root)

from tools.parsers.CSVTools import CSVFile, CSVEntry


class TestCSVFileCreation(unittest.TestCase):

    def test_create_with_legend(self):
        c = CSVFile(legend=["col1", "col2"])
        self.assertIn("col1", c.legend)
        self.assertIn("col2", c.legend)

    def test_create_empty(self):
        c = CSVFile()
        self.assertEqual(len(c.legend), 0)

    def test_duplicate_legend_raises(self):
        with self.assertRaises(ValueError):
            CSVFile(legend=["col1", "col1"])

    def test_custom_separator(self):
        c = CSVFile(separator="\t")
        self.assertEqual(c.separator, "\t")


class TestCSVFileNewLine(unittest.TestCase):

    def test_new_line(self):
        c = CSVFile(legend=["name", "value"])
        line = c.newLine()
        self.assertIsInstance(line, CSVEntry)

    def test_set_and_get_fields(self):
        c = CSVFile(legend=["name", "value"])
        line = c.newLine()
        line["name"] = "test"
        line["value"] = "42"
        self.assertEqual(line["name"], "test")
        self.assertEqual(line["value"], "42")

    def test_case_insensitive_access(self):
        c = CSVFile(legend=["Name", "Value"])
        line = c.newLine()
        line["name"] = "test"
        self.assertEqual(line["Name"], "test")
        self.assertEqual(line["name"], "test")

    def test_nonexistent_column_raises(self):
        c = CSVFile(legend=["name"])
        line = c.newLine()
        with self.assertRaises(KeyError):
            _ = line["nonexistent"]


class TestCSVFileSaveAndParse(unittest.TestCase):

    def setUp(self):
        self.tmpfile = tempfile.mktemp(suffix=".csv")

    def tearDown(self):
        if os.path.exists(self.tmpfile):
            os.remove(self.tmpfile)

    def test_save_and_parse_roundtrip(self):
        c = CSVFile(legend=["col1", "col2"], separator="\t")
        l1 = c.newLine()
        l1["col1"] = "hello"
        l1["col2"] = "world"
        l2 = c.newLine()
        l2["col1"] = "foo"
        l2["col2"] = "bar"
        c.save(self.tmpfile)

        c2 = CSVFile()
        c2.parse(self.tmpfile, separator="\t")
        lines = list(c2)
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[0]["col1"], "hello")
        self.assertEqual(lines[0]["col2"], "world")
        self.assertEqual(lines[1]["col1"], "foo")
        self.assertEqual(lines[1]["col2"], "bar")

    def test_save_comma_separated(self):
        c = CSVFile(legend=["a", "b"], separator=",")
        line = c.newLine()
        line["a"] = "x"
        line["b"] = "y"
        c.save(self.tmpfile)

        c2 = CSVFile()
        c2.parse(self.tmpfile, separator=",")
        self.assertEqual(list(c2)[0]["a"], "x")

    def test_len(self):
        c = CSVFile(legend=["col1"])
        c.newLine()
        c.newLine()
        c.newLine()
        self.assertEqual(len(c), 3)


class TestCSVFileIteration(unittest.TestCase):

    def test_iteration(self):
        c = CSVFile(legend=["name"])
        for i in range(5):
            line = c.newLine()
            line["name"] = str(i)

        names = [line["name"] for line in c]
        self.assertEqual(names, ["0", "1", "2", "3", "4"])

    def test_getitem(self):
        c = CSVFile(legend=["name"])
        line = c.newLine()
        line["name"] = "test"
        self.assertEqual(c[0]["name"], "test")


class TestCSVEntryIteration(unittest.TestCase):

    def test_entry_iteration(self):
        c = CSVFile(legend=["name", "email"])
        line = c.newLine()
        line["name"] = "alice"
        line["email"] = "alice@test.com"
        fields = dict(line)
        self.assertEqual(fields["name"], "alice")
        self.assertEqual(fields["email"], "alice@test.com")


class TestCSVFileAddField(unittest.TestCase):

    def test_add_field(self):
        c = CSVFile(legend=["col1"])
        c.addField("col2")
        self.assertIn("col2", c.legend)

    def test_add_duplicate_field_raises(self):
        c = CSVFile(legend=["col1"])
        with self.assertRaises(ValueError):
            c.addField("col1")


class TestCSVFileStreamToFile(unittest.TestCase):

    def setUp(self):
        self.tmpfile = tempfile.mktemp(suffix=".csv")

    def tearDown(self):
        if os.path.exists(self.tmpfile):
            os.remove(self.tmpfile)

    def test_stream_requires_legend(self):
        c = CSVFile()
        with self.assertRaises(ValueError):
            c.streamToFile(self.tmpfile)

    def test_commit_without_stream_raises(self):
        c = CSVFile(legend=["col1"])
        line = c.newLine()
        with self.assertRaises(ValueError):
            line.commit()


if __name__ == "__main__":
    unittest.main()
