import unittest
import sys
import os

# Add parent of pyGeno to path so we can import tools directly
_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _root not in sys.path:
    sys.path.insert(0, _root)

from tools.UsefulFunctions import (
    complement, reverseComplement, complementTab, reverseComplementTab,
    translateDNA, translateDNA_6Frames,
    findAll,
    encodePolymorphicNucleotide, decodePolymorphicNucleotide,
    decodePolymorphicNucleotide_str,
    getSequenceCombinaisons, polymorphicCodonCombinaisons,
    getNucleotideCodon, showDifferences,
    nucleotides, polymorphicNucleotides, codonTable, translTable,
    UnknownNucleotide,
)


class TestComplement(unittest.TestCase):

    def test_basic_complement(self):
        self.assertEqual(complement("ATCG"), "TAGC")

    def test_complement_is_not_reversed(self):
        self.assertEqual(complement("AAAT"), "TTTA")

    def test_complement_lowercase(self):
        self.assertEqual(complement("atcg"), "tagc")

    def test_complement_mixed_case(self):
        self.assertEqual(complement("AaTt"), "TtAa")

    def test_complement_iupac_codes(self):
        self.assertEqual(complement("R"), "Y")
        self.assertEqual(complement("Y"), "R")
        self.assertEqual(complement("M"), "K")
        self.assertEqual(complement("K"), "M")
        self.assertEqual(complement("W"), "W")
        self.assertEqual(complement("S"), "S")
        self.assertEqual(complement("N"), "N")

    def test_complement_empty(self):
        self.assertEqual(complement(""), "")


class TestReverseComplement(unittest.TestCase):

    def test_basic_rc(self):
        self.assertEqual(reverseComplement("ATCG"), "CGAT")

    def test_rc_palindrome(self):
        self.assertEqual(reverseComplement("AATT"), "AATT")

    def test_rc_single_base(self):
        self.assertEqual(reverseComplement("A"), "T")

    def test_rc_empty(self):
        self.assertEqual(reverseComplement(""), "")


class TestComplementTab(unittest.TestCase):

    def test_basic_list(self):
        result = complementTab(["A", "T", "C", "G"])
        self.assertEqual(result, ["T", "A", "G", "C"])

    def test_empty_string_in_list(self):
        result = complementTab(["A", "", "G"])
        self.assertEqual(result, ["T", "", "C"])

    def test_insertion_in_list(self):
        result = complementTab(["ACT"])
        self.assertEqual(result, [reverseComplement("ACT")])

    def test_empty_list(self):
        self.assertEqual(complementTab([]), [])


class TestReverseComplementTab(unittest.TestCase):

    def test_basic(self):
        result = reverseComplementTab(["A", "T", "C", "G"])
        self.assertEqual(result, ["C", "G", "A", "T"])


class TestTranslateDNA(unittest.TestCase):

    def test_start_codon(self):
        self.assertEqual(translateDNA("ATG"), "M")

    def test_stop_codon_TAA(self):
        self.assertEqual(translateDNA("TAA"), "*")

    def test_stop_codon_TAG(self):
        self.assertEqual(translateDNA("TAG"), "*")

    def test_stop_codon_TGA(self):
        self.assertEqual(translateDNA("TGA"), "*")

    def test_multiple_codons(self):
        self.assertEqual(translateDNA("ATGTTT"), "MF")

    def test_frame_f2(self):
        # f2 skips first nucleotide
        self.assertEqual(translateDNA("AATG", frame="f2"), "M")

    def test_frame_f3(self):
        # f3 skips first two nucleotides
        self.assertEqual(translateDNA("AAATG", frame="f3"), "M")

    def test_frame_r1(self):
        # reverse complement then translate
        seq = "CAT"  # RC = ATG -> M
        self.assertEqual(translateDNA(seq, frame="r1"), "M")

    def test_incomplete_codon_ignored(self):
        # trailing 1-2 nucleotides are ignored
        self.assertEqual(translateDNA("ATGA"), "M")
        self.assertEqual(translateDNA("ATGAT"), "M")

    def test_unknown_frame_raises(self):
        with self.assertRaises(ValueError):
            translateDNA("ATG", frame="x9")

    def test_all_64_codons(self):
        for codon, aa in codonTable.items():
            if '!' not in codon:
                self.assertEqual(translateDNA(codon), aa)

    def test_mitochondrial_table(self):
        # In mt table, AGA is a stop codon
        self.assertEqual(translateDNA("AGA", translTable_id="mt"), "*")
        # In mt table, TGA is W instead of stop
        self.assertEqual(translateDNA("TGA", translTable_id="mt"), "W")


class TestTranslateDNA6Frames(unittest.TestCase):

    def test_returns_6_frames(self):
        result = translateDNA_6Frames("ATGATGATG")
        self.assertEqual(len(result), 6)

    def test_f1_matches_single_call(self):
        seq = "ATGATGATG"
        result = translateDNA_6Frames(seq)
        self.assertEqual(result[0], translateDNA(seq, "f1"))


class TestFindAll(unittest.TestCase):

    def test_basic_find(self):
        self.assertEqual(findAll("ATGATGATG", "ATG"), [0, 3, 6])

    def test_no_match(self):
        self.assertEqual(findAll("AAAA", "CC"), [])

    def test_single_match(self):
        self.assertEqual(findAll("AATGA", "ATG"), [1])

    def test_overlapping_not_found(self):
        # findAll doesn't find overlapping matches
        self.assertEqual(findAll("AAA", "AA"), [0])

    def test_empty_haystack(self):
        self.assertEqual(findAll("", "ATG"), [])


class TestEncodePolymorphicNucleotide(unittest.TestCase):

    def test_single_nucleotide(self):
        self.assertEqual(encodePolymorphicNucleotide("A"), "A")

    def test_two_nucleotides_AG(self):
        self.assertEqual(encodePolymorphicNucleotide("AG"), "R")

    def test_two_nucleotides_CT(self):
        self.assertEqual(encodePolymorphicNucleotide("CT"), "Y")

    def test_two_nucleotides_AC(self):
        self.assertEqual(encodePolymorphicNucleotide("AC"), "M")

    def test_two_nucleotides_TG(self):
        self.assertEqual(encodePolymorphicNucleotide("TG"), "K")

    def test_two_nucleotides_AT(self):
        self.assertEqual(encodePolymorphicNucleotide("AT"), "W")

    def test_two_nucleotides_CG(self):
        self.assertEqual(encodePolymorphicNucleotide("CG"), "S")

    def test_three_nucleotides_CGT(self):
        self.assertEqual(encodePolymorphicNucleotide("CGT"), "B")

    def test_three_nucleotides_AGT(self):
        self.assertEqual(encodePolymorphicNucleotide("AGT"), "D")

    def test_three_nucleotides_ACT(self):
        self.assertEqual(encodePolymorphicNucleotide("ACT"), "H")

    def test_three_nucleotides_ACG(self):
        self.assertEqual(encodePolymorphicNucleotide("ACG"), "V")

    def test_four_nucleotides(self):
        self.assertEqual(encodePolymorphicNucleotide("ACGT"), "N")

    def test_slash_separated(self):
        self.assertEqual(encodePolymorphicNucleotide("A/G"), "R")

    def test_list_input(self):
        self.assertEqual(encodePolymorphicNucleotide(["A", "G"]), "R")

    def test_iupac_input_expands(self):
        # R = A/G, if passed R it should return R
        self.assertEqual(encodePolymorphicNucleotide("R"), "R")


class TestDecodePolymorphicNucleotide(unittest.TestCase):

    def test_decode_R(self):
        self.assertEqual(sorted(decodePolymorphicNucleotide("R")), ["A", "G"])

    def test_decode_Y(self):
        self.assertEqual(sorted(decodePolymorphicNucleotide("Y")), ["C", "T"])

    def test_decode_N(self):
        self.assertEqual(sorted(decodePolymorphicNucleotide("N")), ["A", "C", "G", "T"])

    def test_decode_regular_nucleotide(self):
        self.assertEqual(decodePolymorphicNucleotide("A"), "A")

    def test_decode_invalid_raises(self):
        with self.assertRaises(ValueError):
            decodePolymorphicNucleotide("X")

    def test_str_version(self):
        self.assertEqual(decodePolymorphicNucleotide_str("R"), "A/G")


class TestGetSequenceCombinaisons(unittest.TestCase):

    def test_no_polymorphisms(self):
        self.assertEqual(getSequenceCombinaisons("ATG"), ["ATG"])

    def test_single_polymorphism(self):
        result = sorted(getSequenceCombinaisons("ARG"))
        self.assertEqual(result, sorted(["AAG", "AGG"]))

    def test_multiple_polymorphisms(self):
        result = sorted(getSequenceCombinaisons("RY"))
        expected = sorted(["AC", "AT", "GC", "GT"])
        self.assertEqual(result, expected)


class TestPolymorphicCodonCombinaisons(unittest.TestCase):

    def test_no_ambiguity(self):
        self.assertEqual(polymorphicCodonCombinaisons(["A", "T", "G"]), ["ATG"])

    def test_with_ambiguity(self):
        result = sorted(polymorphicCodonCombinaisons(["R", "T", "G"]))
        self.assertEqual(result, sorted(["ATG", "GTG"]))


class TestGetNucleotideCodon(unittest.TestCase):

    def test_first_codon_pos0(self):
        codon, pos = getNucleotideCodon("ATGCCC", 0)
        self.assertEqual(codon, "ATG")
        self.assertEqual(pos, 0)

    def test_first_codon_pos1(self):
        codon, pos = getNucleotideCodon("ATGCCC", 1)
        self.assertEqual(codon, "ATG")
        self.assertEqual(pos, 1)

    def test_first_codon_pos2(self):
        codon, pos = getNucleotideCodon("ATGCCC", 2)
        self.assertEqual(codon, "ATG")
        self.assertEqual(pos, 2)

    def test_second_codon(self):
        codon, pos = getNucleotideCodon("ATGCCC", 3)
        self.assertEqual(codon, "CCC")
        self.assertEqual(pos, 0)

    def test_out_of_range(self):
        self.assertIsNone(getNucleotideCodon("ATG", 5))

    def test_negative_index(self):
        self.assertIsNone(getNucleotideCodon("ATG", -1))


class TestShowDifferences(unittest.TestCase):

    def test_identical_sequences(self):
        result = showDifferences("ATG", "ATG")
        self.assertIn("-", result)
        self.assertNotIn("|", result)

    def test_different_sequences(self):
        result = showDifferences("ATG", "ATC")
        self.assertIn("|", result)

    def test_different_lengths(self):
        result = showDifferences("ATGC", "AT")
        self.assertIn("#", result)


if __name__ == "__main__":
    unittest.main()
