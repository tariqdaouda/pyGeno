from . import configuration as conf

from .pyGenoObjectBases import *

import rabaDB.fields as rf

from .tools import UsefulFunctions as uf
from .Exon import *
from .SNP import SNP_INDEL

from .tools.BinarySequence import NucBinarySequence


class Transcript_Raba(pyGenoRabaObject) :
	"""The wrapped Raba object that really holds the data"""
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	id = rf.Primitive()
	name = rf.Primitive()
	length = rf.Primitive()
	start = rf.Primitive()
	end = rf.Primitive()
	coding = rf.Primitive()
	biotype = rf.Primitive()
	selenocysteine = rf.RList()
	
	genome = rf.RabaObject('Genome_Raba')
	chromosome = rf.RabaObject('Chromosome_Raba')
	gene = rf.RabaObject('Gene_Raba')
	protein = rf.RabaObject('Protein_Raba')
	exons = rf.Relation('Exon_Raba')
	
	def _curate(self) :
		if self.name != None :
			self.name = self.name.upper()
		
		self.length = abs(self.end - self.start)
		have_CDS_start = False
		have_CDS_end = False
		for exon in self.exons :
			if exon.CDS_start is not None :
				have_CDS_start = True
			if exon.CDS_end is not None :
				have_CDS_end = True

		if have_CDS_start and have_CDS_end :
			self.coding = True
		else :
			self.coding = False

class Transcript(pyGenoRabaObjectWrapper) :
	"""The wrapper for playing with Transcripts"""
	
	_wrapped_class = Transcript_Raba

	def __init__(self, *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)
		self.exons = RLWrapper(self, Exon, self.wrapped_object.exons)
		self._load_sequencesTriggers = set(["UTR5", "UTR3", "cDNA", "sequence", "data"])
		self.exonsDict = {}
	
	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		if issubclass(objectType, SNP_INDEL) :
			# conf.db.enableDebug(True)
			f = RabaQuery(objectType, namespace = self._wrapped_class._raba_namespace)
			coolArgs['species'] = self.genome.species
			coolArgs['chromosomeNumber'] = self.chromosome.number
			coolArgs['start >='] = self.start
			coolArgs['start <'] = self.end
		
			if len(args) > 0 and type(args[0]) is list :
				for a in args[0] :
					if type(a) is dict :
						f.addFilter(**a)
			else :
				f.addFilter(*args, **coolArgs)

			return f
		
		return pyGenoRabaObjectWrapper._makeLoadQuery(self, objectType, *args, **coolArgs)
		
	def _load_data(self) :
		def getV(k) :
			return pyGenoRabaObjectWrapper.__getattribute__(self, k)

		def setV(k, v) :
			return pyGenoRabaObjectWrapper.__setattr__(self, k, v)

		self.data = []
		cDNA = []
		UTR5 = []
		UTR3 = []
		exons = []
		prime5 = True
		for ee in self.wrapped_object.exons :
			e = pyGenoRabaObjectWrapper_metaclass._wrappers[Exon_Raba](wrapped_object_and_bag = (ee, getV('bagKey')))
			self.exonsDict[(e.start, e.end)] = e
			exons.append(e)
			self.data.extend(e.data)
			
			if e.hasCDS() :
				UTR5.append(''.join(e.UTR5))

				e._patch_seleno(e, self.selenocysteine)

				if len(cDNA) == 0 and e.frame != 0 :
					e.CDS = e.CDS[e.frame:]
					
					if e.strand == '+':
						e.CDS_start += e.frame
					else:
						e.CDS_end -= e.frame
				
				if len(e.CDS):
					cDNA.append(''.join(e.CDS))
				UTR3.append(''.join(e.UTR3))
				prime5 = False
			else :
				if prime5 :
					UTR5.append(''.join(e.data))
				else :
					UTR3.append(''.join(e.data))
		
		sequence = ''.join(self.data)
		cDNA = ''.join(cDNA)
		UTR5 = ''.join(UTR5)
		UTR3 = ''.join(UTR3)
		setV('exons', exons)
		setV('sequence', sequence)
		setV('cDNA', cDNA)
		setV('UTR5', UTR5)
		setV('UTR3', UTR3)
		
		if len(cDNA) > 0 and len(cDNA) % 3 != 0 :
			setV('flags', {'DUBIOUS' : True, 'cDNA_LEN_NOT_MULT_3': True})
		else :
			setV('flags', {'DUBIOUS' : False, 'cDNA_LEN_NOT_MULT_3': False})

	def _load_bin_sequence(self) :
		self.bin_sequence = NucBinarySequence(self.sequence)
		self.bin_UTR5 =  NucBinarySequence(self.UTR5)
		self.bin_cDNA =  NucBinarySequence(self.cDNA.replace('!', 'T' ))
		self.bin_UTR3 =  NucBinarySequence(self.UTR3)

	def getNucleotideCodon(self, cdnaX1) :
		"""Returns the entire codon of the nucleotide at pos cdnaX1 in the cdna, and the position of that nocleotide in the codon"""
		return uf.getNucleotideCodon(self.cDNA, cdnaX1)

	def getCodon(self, i) :
		"""returns the ith codon"""
		return self.getNucleotideCodon(i*3)[0]

	def iterCodons(self) :
		"""iterates through the codons"""
		for i in range(len(self.cDNA)/3) :
			yield self.getCodon(i)

	def find(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.bin_sequence.find(sequence)

	def findAll(self, sequence):
		"""Returns a list of all positions where sequence was found"""
		return self.bin_sequence.findAll(sequence)

	def findIncDNA(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.bin_cDNA.find(sequence)

	def findAllIncDNA(self, sequence) :
		"""Returns a list of all positions where sequence was found in the cDNA"""
		return self.bin_cDNA.findAll(sequence)

	def getcDNALength(self) :
		"""returns the length of the cDNA"""
		return len(self.cDNA)

	def findInUTR5(self, sequence) :
		"""return the position of the first occurance of sequence in the 5'UTR"""
		return self.bin_UTR5.find(sequence)

	def findAllInUTR5(self, sequence) :
		"""Returns a list of all positions where sequence was found in the 5'UTR"""
		return self.bin_UTR5.findAll(sequence)

	def getUTR5Length(self) :
		"""returns the length of the 5'UTR"""
		return len(self.bin_UTR5)

	def findInUTR3(self, sequence) :
		"""return the position of the first occurance of sequence in the 3'UTR"""
		return self.bin_UTR3.find(sequence)

	def findAllInUTR3(self, sequence) :
		"""Returns a lits of all positions where sequence was found in the 3'UTR"""
		return self.bin_UTR3.findAll(sequence)

	def getUTR3Length(self) :
		"""returns the length of the 3'UTR"""
		return len(self.bin_UTR3)

	def getNbCodons(self) :
		"""returns the number of codons in the transcript"""
		return len(self.cDNA)/3
	
	def __getattribute__(self, name) :
		return pyGenoRabaObjectWrapper.__getattribute__(self, name)

	def __getitem__(self, i) :
		return self.sequence[i]

	def __len__(self) :
		return len(self.sequence)

	def __str__(self) :
		return """Transcript, id: %s, name: %s > %s""" %(self.id, self.name, str(self.gene))
