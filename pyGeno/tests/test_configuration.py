import unittest
import sys
import os
from unittest.mock import MagicMock
from configparser import ConfigParser

# Mock rabaDB before importing configuration
if 'rabaDB' not in sys.modules:
    sys.modules['rabaDB'] = MagicMock()
    sys.modules['rabaDB.rabaSetup'] = MagicMock()
    sys.modules['rabaDB.Raba'] = MagicMock()

from pyGeno.configuration import checkPythonVersion, version, prettyVersion


class TestConfigParserCompat(unittest.TestCase):
    """Verify that ConfigParser (not SafeConfigParser) works correctly."""

    def test_configparser_read_write(self):
        parser = ConfigParser()
        parser.add_section("test_section")
        parser.set("test_section", "key", "value")
        self.assertEqual(parser.get("test_section", "key"), "value")

    def test_configparser_items(self):
        parser = ConfigParser()
        parser.add_section("section")
        parser.set("section", "a", "1")
        parser.set("section", "b", "2")
        items = dict(parser.items("section"))
        self.assertEqual(items["a"], "1")
        self.assertEqual(items["b"], "2")


class TestPythonVersionCheck(unittest.TestCase):

    def test_python_version_is_valid(self):
        result = checkPythonVersion()
        self.assertTrue(result)


class TestVersionFunctions(unittest.TestCase):

    def test_version_returns_tuple(self):
        v = version()
        self.assertIsInstance(v, tuple)
        self.assertEqual(len(v), 6)

    def test_pretty_version_returns_string(self):
        pv = prettyVersion()
        self.assertIsInstance(pv, str)
        self.assertIn("pyGeno", pv)


if __name__ == "__main__":
    unittest.main()
