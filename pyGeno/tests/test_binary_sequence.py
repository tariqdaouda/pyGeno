import unittest
import sys
import os

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _root not in sys.path:
    sys.path.insert(0, _root)

from tools.BinarySequence import NucBinarySequence, AABinarySequence


class TestNucBinarySequenceEncode(unittest.TestCase):

    def test_simple_sequence(self):
        bs = NucBinarySequence("ATCG")
        self.assertEqual(len(bs), 4)
        self.assertEqual(bs.getDefaultSequence(), "ATCG")

    def test_empty_polymorphisms(self):
        bs = NucBinarySequence("ATCG")
        self.assertEqual(bs.getPolymorphisms(), [])

    def test_polymorphic_sequence(self):
        bs = NucBinarySequence("R")  # R = A/G
        polys = bs.getPolymorphisms()
        self.assertEqual(len(polys), 1)
        self.assertEqual(polys[0][0], 0)  # position 0
        self.assertIn("A", polys[0][1])
        self.assertIn("G", polys[0][1])

    def test_polymorphic_n(self):
        bs = NucBinarySequence("N")  # N = A/C/G/T
        polys = bs.getPolymorphisms()
        self.assertEqual(len(polys), 1)
        self.assertEqual(len(polys[0][1]), 4)

    def test_mixed_sequence(self):
        bs = NucBinarySequence("ATRGC")
        self.assertEqual(len(bs), 5)
        polys = bs.getPolymorphisms()
        self.assertEqual(len(polys), 1)
        self.assertEqual(polys[0][0], 2)  # R is at position 2

    def test_default_sequence_with_polymorphism(self):
        bs = NucBinarySequence("ARC")
        default = bs.getDefaultSequence()
        self.assertEqual(len(default), 3)
        # Default takes last allele of polymorphism
        self.assertIn(default[1], "AG")


class TestNucBinarySequenceFind(unittest.TestCase):

    def test_find_exact(self):
        bs = NucBinarySequence("ATGATGATG")
        pos = bs.find("ATG")
        self.assertEqual(pos, 0)

    def test_find_middle(self):
        bs = NucBinarySequence("CCATGCC")
        pos = bs.find("ATG")
        self.assertEqual(pos, 2)

    def test_find_not_found(self):
        bs = NucBinarySequence("AAAA")
        pos = bs.find("GGG")
        self.assertEqual(pos, -1)

    def test_find_all(self):
        bs = NucBinarySequence("ATGATGATG")
        positions = bs.findAll("ATG")
        self.assertEqual(positions, [0, 3, 6])

    def test_find_all_no_match(self):
        bs = NucBinarySequence("AAAA")
        positions = bs.findAll("GGG")
        self.assertEqual(positions, [])


class TestNucBinarySequenceVariants(unittest.TestCase):

    def test_no_polymorphism_variants(self):
        bs = NucBinarySequence("ATCG")
        stopped, variants = bs.getSequenceVariants()
        self.assertFalse(stopped)
        self.assertEqual(variants, ["ATCG"])

    def test_single_polymorphism_variants(self):
        bs = NucBinarySequence("ARC")  # R = A/G
        stopped, variants = bs.getSequenceVariants()
        self.assertFalse(stopped)
        self.assertEqual(sorted(variants), sorted(["AAC", "AGC"]))

    def test_n_polymorphism_variants(self):
        bs = NucBinarySequence("N")  # 4 variants
        stopped, variants = bs.getSequenceVariants()
        self.assertFalse(stopped)
        self.assertEqual(len(variants), 4)

    def test_max_variant_limit(self):
        # NN = 16 variants, set limit to 4
        bs = NucBinarySequence("NN")
        stopped, variants = bs.getSequenceVariants(maxVariantNumber=4)
        self.assertTrue(stopped)

    def test_nb_variants(self):
        bs = NucBinarySequence("ARN")  # R=2 * N=4 = 8
        self.assertEqual(bs.getNbVariants(0), 8)

    def test_nb_variants_no_poly(self):
        bs = NucBinarySequence("ATCG")
        self.assertEqual(bs.getNbVariants(0), 1)


class TestNucBinarySequenceLen(unittest.TestCase):

    def test_length(self):
        bs = NucBinarySequence("ATCG")
        self.assertEqual(len(bs), 4)

    def test_length_with_polymorphism(self):
        bs = NucBinarySequence("ARCG")  # R expands but sequence is still 4
        self.assertEqual(len(bs), 4)


class TestNucBinarySequenceFindPolymorphisms(unittest.TestCase):

    def test_identical_no_polymorphisms(self):
        bs = NucBinarySequence("ATCG")
        result = bs.findPolymorphisms("ATCG")
        self.assertEqual(result, [])

    def test_single_mismatch(self):
        bs = NucBinarySequence("ATCG")
        result = bs.findPolymorphisms("AACG")
        self.assertEqual(result, [1])


class TestAABinarySequence(unittest.TestCase):

    def test_simple_sequence(self):
        bs = AABinarySequence("MLPK")
        self.assertEqual(len(bs), 4)
        self.assertEqual(bs.getDefaultSequence(), "MLPK")

    def test_find_in_protein(self):
        bs = AABinarySequence("MLPKADEW")
        pos = bs.find("ADE")
        self.assertEqual(pos, 4)

    def test_polymorphic_protein(self):
        bs = AABinarySequence("M/LPK")
        polys = bs.getPolymorphisms()
        self.assertEqual(len(polys), 1)
        self.assertIn("M", polys[0][1])
        self.assertIn("L", polys[0][1])


class TestBinarySequenceGetItem(unittest.TestCase):

    def test_getitem(self):
        bs = NucBinarySequence("ATCG")
        # Each element should be an integer
        for i in range(4):
            self.assertIsInstance(bs[i], int)

    def test_setitem(self):
        bs = NucBinarySequence("ATCG")
        original = bs[0]
        bs[0] = 0
        self.assertEqual(bs[0], 0)
        bs[0] = original

    def test_getchar(self):
        bs = NucBinarySequence("ATCG")
        self.assertEqual(bs.getChar(0), "A")
        self.assertEqual(bs.getChar(1), "T")
        self.assertEqual(bs.getChar(2), "C")
        self.assertEqual(bs.getChar(3), "G")


if __name__ == "__main__":
    unittest.main()
