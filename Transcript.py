import configuration as conf

from pyGenoObjectBases import *

import rabaDB.fields as rf

from tools import UsefulFunctions as uf
#from Protein import Protein
from Exon import *

from tools.BinarySequence import NucBinarySequence


class Transcript_Raba(pyGenoRabaObject) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	id = rf.Primitive()
	name = rf.Primitive()
	length = rf.Primitive()
	start = rf.Primitive()
	end = rf.Primitive()
	coding = rf.Primitive()
	
	genome = rf.RabaObject('Genome_Raba')
	chromosome = rf.RabaObject('Chromosome_Raba')
	gene = rf.RabaObject('Gene_Raba')
	protein = rf.RabaObject('Protein_Raba')
	exons = rf.Relation('Exon_Raba')
	
	def _curate(self) :
		if self.name != None :
			self.name = self.name.upper()
		
		self.start = self.exons[0].start
		self.end = self.exons[-1].end
		self.length = abs(self.end - self.start)
		try :
			self.CDS_start = self.exons[0].CDS_start
			self.CDS_end = self.exons[-1].CDS_end
			self.CDS_length = abs(self.CDS_end - self.CDS_start)
			self.coding = True
		except :
			self.CDS_start, self.CDS_end = None, None
			self.CDS_length = None
			self.coding = False

class Transcript(pyGenoRabaObjectWrapper) :

	_wrapped_class = Transcript_Raba

	def __init__(self, *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)
		self.loadedSequences = False
		self.loadBinarySequence = True

	def _loadSequences(self) :
		if not pyGenoRabaObjectWrapper.__getattribute__(self, 'loadedSequences') :
			def getV(k) :
				return pyGenoRabaObjectWrapper.__getattribute__(self, k)

			def setV(k, v) :
				return pyGenoRabaObjectWrapper.__setattr__(self, k, v)

			sequence = []
			CDNA = []
			UTR5 = []
			UTR3 = []
			exons = []
			prime5 = True
			for ee in self.wrapped_object.exons :
				e = pyGenoRabaObjectWrapper_metaclass._wrappers[Exon_Raba](wrapped_object_and_bag = (ee, getV('bagKey')))
				exons.append(e)
				sequence.append(e.sequence)

				if e.hasCDS() :
					UTR5.append(e.UTR5)
					CDNA.append(e.CDS)
					UTR3.append(e.UTR3)
					prime5 = False
				else :
					if prime5 :
						UTR5.append(e.sequence)
					else :
						UTR3.append(e.sequence)

			setV('exons', exons)
			setV('sequence', ''.join(sequence))
			setV('CDNA', ''.join(CDNA))
			setV('UTR5', ''.join(UTR5))
			setV('UTR3', ''.join(UTR3))

			if len(getV('CDNA')) % 3 != 0 :
				setV('flags', {'DUBIOUS' : True, 'CDNA_LEN_NOT_MULT_3': True})
			else :
				setV('flags', {'DUBIOUS' : False, 'CDNA_LEN_NOT_MULT_3': False})

			setV('loadedSequences', True)

	def _load_bin_sequence(self) :
		self.bin_sequence = NucBinarySequence(self.sequence)
		self.bin_UTR5 =  NucBinarySequence(self.UTR5)
		self.bin_CDS =  NucBinarySequence(self.CDS)
		self.bin_UTR3 =  NucBinarySequence(self.UTR3)

	def getNucleotideCodon(self, cdnaX1) :
		"Returns the entire codon of the nucleotide at pos cdnaX1 in the cdna, and the position of that nocleotide in the codon"
		return uf.getNucleotideCodon(self.CDNA, cdnaX1)

	def getCodon(self, i) :
		"returns the ith codon"
		return self.getNucleotideCodon(i*3)[0]

	def iterCodons(self) :
		"iterates through the codons"
		for i in range(len(self.CDNA)/3) :
			yield self.getCodon(i)

	def find(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.bin_Sequence.find(sequence)

	def findAll(self, seqence):
		"""Returns a lits of all positions where sequence was found"""
		return self.bin_Sequence.findAll(sequence)

	def findInCDNA(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.bin_CDNA.find(sequence)

	def findAllInCDNA(self, seqence):
		"""Returns a lits of all positions where sequence was found"""
		return self.bin_CDNA.findAll(sequence)

	def getCDNALength(self):
		return len(self.CDNA)

	def findInUTR5(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.bin_UTR5.find(sequence)

	def findAllInUTR5(self, seqence):
		"""Returns a lits of all positions where sequence was found"""
		return self.bin_UTR5.findAll(sequence)

	def getUTR5Length(self):
		return len(self.bin_UTR5)

	def findInUTR3(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.bin_UTR3.find(sequence)

	def findAllInUTR3(self, seqence):
		"""Returns a lits of all positions where sequence was found"""
		return self.bin_UTR3.findAll(sequence)

	def getUTR3Length(self):
		return len(self.bin_UTR3)

	#<7iyed>
	def getCodonAffinityMap(self, chunkRatio = 0.05) :
		chunks = []
		if len(self.CDNA) < 3 :
			return None

		for i in range(int(1/chunkRatio)) :
			chunks.append({'low_aff':0, 'high_aff':0})

		for i in range(self.getNbCodons()) :
			c = int(N.floor(i/(self.getNbCodons()*chunkRatio)))
			codon = self.getCodon(i)[0]
			if codon in uf.codonTable.keys():
				if uf.codonAffinity[codon] == 'low' :
					chunks[c]['low_aff'] += 1
				elif uf.codonAffinity[codon] == 'high' :
					chunks[c]['high_aff'] += 1
			else :
				for polyCodon in uf.polymorphicCodonCombinaisons(codon) :
					if uf.codonAffinity[polyCodon] == 'low' :
						chunks[c]['low_aff'] += 1
					elif uf.codonAffinity[polyCodon] == 'high' :
						chunks[c]['high_aff'] += 1

		return chunks

	def getCodonUsage(self) :
		codonUsage = {}
		for k in uf.codonTable.keys() :
			codonUsage[k] = 0

		for i in range(self.getNbCodons()) :
			codon = self.getCodon(i)[0]
			if codon in uf.codonTable.keys():
				codonUsage[codon] += 1
			else :
				for polyCodon in uf.polymorphicCodonCombinaisons(codon) :
					codonUsage[polyCodon] += 1
		return codonUsage
	#</7iyed>

	def pluck(self):
		"""Returns a plucked object. Plucks the transcript off the tree, set the value of self.gene into str(self.gene). This effectively disconnects the object and
		makes it much more lighter in case you'd like to pickle it"""
		e = copy.copy(self)
		e.gene = str(self.gene)
		return e

	def getNbCodons(self) :
		return len(self.CDNA)/3

	def __getattribute__(self, name) :
		try :
			return pyGenoRabaObjectWrapper.__getattribute__(self, name)
		except AttributeError :
			pyGenoRabaObjectWrapper.__getattribute__(self, '_loadSequences')()

		return pyGenoRabaObjectWrapper.__getattribute__(self, name)

	def __getitem__(self, i) :
		return self.sequence[i]

	def __len__(self) :
		return len(self.sequence)

	def __str__(self) :
		return """Transcript, id: %s, name: %s > %s""" %(self.id, self.name, str(self.gene))

