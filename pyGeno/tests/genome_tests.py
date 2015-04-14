import unittest
from pyGeno.Genome import *

import pyGeno.bootstrap as B


class pyGenoSNPTests(unittest.TestCase):

	def setUp(self):
		self.ref = Genome(name = 'GRCh37.75_Y-Only')

	def tearDown(self):
		pass

	def test_vanilla(self) :
		dummy = Genome(name = 'GRCh37.75_Y-Only', SNPs = 'dummySRY_AGN')
		persProt = dummy.get(Protein, id = 'ENSP00000438917')[0]
		refProt = self.ref.get(Protein, id = 'ENSP00000438917')[0]
		
		self.assertEqual('ATGCAATCATATGCTTCTGC', refProt.transcript.cDNA[:20])
		self.assertEqual('HTGCAATCATATGCTTCTGC', persProt.transcript.cDNA[:20])
		
	def test_noModif(self) :
		from pyGeno.SNPFiltering import SNPFilter

		class MyFilter(SNPFilter) :
			def __init__(self) :
				SNPFilter.__init__(self)

			def filter(self, chromosome, dummySRY_AGN) :
				return None

		dummy = Genome(name = 'GRCh37.75_Y-Only', SNPs = 'dummySRY_AGN', SNPFilter = MyFilter())
		persProt = dummy.get(Protein, id = 'ENSP00000438917')[0]
		refProt = self.ref.get(Protein, id = 'ENSP00000438917')[0]

		self.assertEqual(persProt.transcript.cDNA[:20], refProt.transcript.cDNA[:20])
	
	def test_insert(self) :
		from pyGeno.SNPFiltering import SNPFilter

		class MyFilter(SNPFilter) :
			def __init__(self) :
				SNPFilter.__init__(self)

			def filter(self, chromosome, dummySRY_AGN) :
				from pyGeno.SNPFiltering import  SequenceInsert
					
				refAllele = chromosome.refSequence[dummySRY_AGN.start]
				return SequenceInsert('TCA')

		dummy = Genome(name = 'GRCh37.75_Y-Only', SNPs = 'dummySRY_AGN', SNPFilter = MyFilter())
		persProt = dummy.get(Protein, id = 'ENSP00000438917')[0]
		refProt = self.ref.get(Protein, id = 'ENSP00000438917')[0]

		self.assertEqual('ATGCAATCATATGCTTCTGC', refProt.transcript.cDNA[:20])
		self.assertEqual('ATGATGCAATCATATGCTTC', persProt.transcript.cDNA[:20])
	
	def test_SNP(self) :
		from pyGeno.SNPFiltering import SNPFilter

		class MyFilter(SNPFilter) :
			def __init__(self) :
				SNPFilter.__init__(self)

			def filter(self, chromosome, dummySRY_AGN) :
				from pyGeno.SNPFiltering import SequenceSNP
				
				return SequenceSNP(dummySRY_AGN.alt)
		
		dummy = Genome(name = 'GRCh37.75_Y-Only', SNPs = 'dummySRY_AGN', SNPFilter = MyFilter())
		persProt = dummy.get(Protein, id = 'ENSP00000438917')[0]
		refProt = self.ref.get(Protein, id = 'ENSP00000438917')[0]
		self.assertEqual('M', refProt.sequence[0])
		self.assertEqual('L', persProt.sequence[0])
		
	def test_deletion(self) :
		from pyGeno.SNPFiltering import SNPFilter

		class MyFilter(SNPFilter) :
			def __init__(self) :
				SNPFilter.__init__(self)

			def filter(self, chromosome, dummySRY_AGN) :
				from pyGeno.SNPFiltering import SequenceDel
					
				refAllele = chromosome.refSequence[dummySRY_AGN.start]
				return SequenceDel(1)
		
		dummy = Genome(name = 'GRCh37.75_Y-Only', SNPs = 'dummySRY_AGN', SNPFilter = MyFilter())
		persProt = dummy.get(Protein, id = 'ENSP00000438917')[0]
		refProt = self.ref.get(Protein, id = 'ENSP00000438917')[0]

		self.assertEqual('ATGCAATCATATGCTTCTGC', refProt.transcript.cDNA[:20])
		self.assertEqual('TGCAATCATATGCTTCTGCT', persProt.transcript.cDNA[:20])

	def test_bags(self) :
		dummy = Genome(name = 'GRCh37.75_Y-Only')
		self.assertEqual(dummy.wrapped_object, self.ref.wrapped_object)
	
	def test_prot_find(self) :
		prot = self.ref.get(Protein, id = 'ENSP00000438917')[0]
		needle = prot.sequence[:10]
		self.assertEqual(0, prot.find(needle))
		needle = prot.sequence[-10:]
		self.assertEqual(len(prot)-10, prot.find(needle))

	def test_trans_find(self) :
		trans = self.ref.get(Transcript, name = "SRY-001")[0]
		self.assertEqual(0, trans.find(trans[:5]))
	
	def test_import_remote_genome(self) :
		self.assertRaises(KeyError, B.importRemoteGenome, "Human.GRCh37.75_Y-Only.tar.gz")

	def test_import_remote_snps(self) :
		self.assertRaises(KeyError, B.importRemoteSNPs, "Human_agnostic.dummySRY.tar.gz")

def runTests() :
	try :
		B.importGenome("Human.GRCh37.75_Y-Only.tar.gz")
	except KeyError :
		print "--> Seems to already exist in db"
		
	try :
		B.importSNPs("Human_agnostic.dummySRY.tar.gz")
	except KeyError :
		print "--> Seems to already exist in db"
		
	unittest.main()

if __name__ == "__main__" :
	runTests()
