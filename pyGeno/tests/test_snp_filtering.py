import unittest
import sys
import os
from unittest.mock import MagicMock

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _root not in sys.path:
    sys.path.insert(0, _root)

# Mock rabaDB before importing SNPFiltering (which imports configuration -> rabaDB)
if 'rabaDB' not in sys.modules:
    sys.modules['rabaDB'] = MagicMock()
    sys.modules['rabaDB.rabaSetup'] = MagicMock()
    sys.modules['rabaDB.Raba'] = MagicMock()

from pyGeno.SNPFiltering import (
    SequenceSNP, SequenceInsert, SequenceDel,
    SNPFilter, DefaultSNPFilter,
    Sequence_modifiers,
)


class TestSequenceModifiers(unittest.TestCase):

    def test_init_default_sources(self):
        sm = Sequence_modifiers()
        self.assertIsInstance(sm.sources, dict)

    def test_add_source(self):
        sm = Sequence_modifiers(sources={})
        sm.addSource("snp1", "some_data")
        self.assertEqual(sm.sources["snp1"], "some_data")


class TestSequenceSNP(unittest.TestCase):

    def test_single_allele(self):
        snp = SequenceSNP("A")
        self.assertEqual(snp.alleles, "A")

    def test_multiple_alleles_string(self):
        snp = SequenceSNP("AG")
        self.assertEqual(snp.alleles, "R")

    def test_multiple_alleles_list(self):
        snp = SequenceSNP(["A", "G"])
        self.assertEqual(snp.alleles, "R")

    def test_ct_alleles(self):
        snp = SequenceSNP("CT")
        self.assertEqual(snp.alleles, "Y")

    def test_all_four_alleles(self):
        snp = SequenceSNP("ACGT")
        self.assertEqual(snp.alleles, "N")

    def test_sources_passed(self):
        snp = SequenceSNP("A", sources={"src": "val"})
        self.assertEqual(snp.sources["src"], "val")


class TestSequenceInsert(unittest.TestCase):

    def test_basic_insert(self):
        ins = SequenceInsert("ACTG")
        self.assertEqual(ins.bases, "ACTG")
        self.assertEqual(ins.offset, 0)

    def test_insert_with_ref_prefix(self):
        # Format like C/CCTGGAA (dbSNP style)
        ins = SequenceInsert("CCTGGAA", ref="C")
        self.assertEqual(ins.bases, "CTGGAA")
        self.assertEqual(ins.offset, 0)  # offset = len(ref) - 1 = 0

    def test_insert_with_longer_ref_prefix(self):
        # Format like CCT/CCTGGAA (samtools style)
        ins = SequenceInsert("CCTGGAA", ref="CCT")
        self.assertEqual(ins.bases, "GGAA")
        self.assertEqual(ins.offset, 2)  # offset = len(ref) - 1 = 2

    def test_insert_bad_ref_raises(self):
        with self.assertRaises(NotImplementedError):
            SequenceInsert("CCTGGAA", ref="GGG")

    def test_insert_default_ref(self):
        ins = SequenceInsert("XXX")
        self.assertEqual(ins.bases, "XXX")
        self.assertEqual(ins.offset, 0)


class TestSequenceDel(unittest.TestCase):

    def test_basic_deletion(self):
        d = SequenceDel(5)
        self.assertEqual(d.length, 5)
        self.assertEqual(d.offset, 0)

    def test_deletion_with_alt_prefix(self):
        # Format like CCTGGAA/C (dbSNP style)
        d = SequenceDel(7, ref="CCTGGAA", alt="C")
        self.assertEqual(d.offset, 1)
        self.assertEqual(d.length, 6)

    def test_deletion_bad_alt_raises(self):
        with self.assertRaises(NotImplementedError):
            SequenceDel(7, ref="CCTGGAA", alt="GGG")

    def test_deletion_alt_without_ref_raises(self):
        with self.assertRaises(Exception):
            SequenceDel(5, alt="C")

    def test_deletion_default_alt(self):
        d = SequenceDel(3)
        self.assertEqual(d.length, 3)
        self.assertEqual(d.offset, 0)


class TestSNPFilter(unittest.TestCase):

    def test_abstract_filter_raises(self):
        f = SNPFilter()
        with self.assertRaises(NotImplementedError):
            f.filter(None)

    def test_default_filter_instantiation(self):
        f = DefaultSNPFilter()
        self.assertIsInstance(f, SNPFilter)


class TestPython314Compatibility(unittest.TestCase):
    """Tests specifically verifying Python 3.14 compatibility fixes."""

    def test_not_implemented_error_is_raised_not_not_implemented(self):
        """Verify NotImplementedError (exception) is raised, not NotImplemented (constant)."""
        with self.assertRaises(NotImplementedError):
            SequenceInsert("CCTGGAA", ref="GGG")

        with self.assertRaises(NotImplementedError):
            SequenceDel(7, ref="CCTGGAA", alt="GGG")

        f = SNPFilter()
        with self.assertRaises(NotImplementedError):
            f.filter(None)


if __name__ == "__main__":
    unittest.main()
