from .pyGenoObjectBases import *
from .SNP import SNP_INDEL

import rabaDB.fields as rf
from .tools import UsefulFunctions as uf
from .tools.BinarySequence import NucBinarySequence

class Exon_Raba(pyGenoRabaObject) :
	"""The wrapped Raba object that really holds the data"""
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	id = rf.Primitive()
	number = rf.Primitive()
	start = rf.Primitive()
	end = rf.Primitive()
	length = rf.Primitive()
	CDS_length = rf.Primitive()
	CDS_start = rf.Primitive()
	CDS_end = rf.Primitive()
	frame = rf.Primitive()
	strand = rf.Primitive()

	genome = rf.RabaObject('Genome_Raba')
	chromosome = rf.RabaObject('Chromosome_Raba')
	gene = rf.RabaObject('Gene_Raba')
	transcript = rf.RabaObject('Transcript_Raba')
	protein = rf.RabaObject('Protein_Raba')

	def _curate(self) :
		if self.start != None and self.end != None :
			if self.start > self.end :
				self.start, self.end = self.end, self.start
			self.length = self.end-self.start

		if self.CDS_start != None and self.CDS_end != None :
			if self.CDS_start > self.CDS_end :
				self.CDS_start, self.CDS_end = self.CDS_end, self.CDS_start
			self.CDS_length = self.CDS_end - self.CDS_start
		
		if self.number != None :
			self.number = int(self.number)

		if not self.frame or self.frame == '.' :
			self.frame = None
		else :
			self.frame = int(self.frame)

class Exon(pyGenoRabaObjectWrapper) :
	"""The wrapper for playing with Exons"""
		
	_wrapped_class = Exon_Raba

	def __init__(self, *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)
		self._load_sequencesTriggers = set(["UTR5", "UTR3", "CDS", "sequence", "data"])

	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		if issubclass(objectType, SNP_INDEL) :
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
		data = self.chromosome.getSequenceData(slice(self.start,self.end))

		diffLen = (self.end-self.start) - len(data)
		
		if self.strand == '+' :
			self.data = data
		else :
			self.data = uf.reverseComplementTab(data)

		if self.hasCDS() :
			start = self.CDS_start-self.start
			end = self.CDS_end-self.start
			
			if self.strand == '+' :
				self.UTR5 = self.data[:start]
				self.CDS = self.data[start:end+diffLen]
				self.UTR3 = self.data[end+diffLen:]
			else :
				self.UTR5 = self.data[:len(self.data)-(end-diffLen)]
				self.CDS = self.data[len(self.data)-(end-diffLen):len(self.data)-start]
				self.UTR3 = self.data[len(self.data)-start:]
		else :
			self.UTR5 = ''
			self.CDS = ''
			self.UTR3 = ''

		self.sequence = ''.join(self.data)

	def _load_bin_sequence(self) :
		self.bin_sequence = NucBinarySequence(self.sequence)
		self.bin_UTR5 =  NucBinarySequence(self.UTR5)
		self.bin_CDS =  NucBinarySequence(self.CDS)
		self.bin_UTR3 =  NucBinarySequence(self.UTR3)
		
	def _patch_seleno(self, e, selenocysteine):
		if selenocysteine is not None:
			for position in selenocysteine:
				if e.CDS_start <= position <= e.CDS_end:

					if e.strand == '+':
						ajusted_position = position - e.CDS_start
					else:
						ajusted_position = e.CDS_end - position - 3

					if e.CDS[ajusted_position] == 'T':
						e.CDS = list(e.CDS)
						e.CDS[ajusted_position] = '!'

	def hasCDS(self) :
		"""returns true or false depending on if the exon has a CDS"""
		if self.CDS_start != None and self.CDS_end != None:
			return True
		return False

	def getCDSLength(self) :
		"""returns the length of the CDS sequence"""
		return len(self.CDS)

	def find(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.bin_sequence.find(sequence)

	def findAll(self, sequence):
		"""Returns a lits of all positions where sequence was found"""
		return self.bin_sequence.findAll(sequence)

	def findInCDS(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.bin_CDS.find(sequence)

	def findAllInCDS(self, sequence):
		"""Returns a lits of all positions where sequence was found"""
		return self.bin_CDS.findAll(sequence)

	def pluck(self) :
		"""Returns a plucked object. Plucks the exon off the tree, set the value of self.transcript into str(self.transcript). This effectively disconnects the object and
		makes it much more lighter in case you'd like to pickle it"""
		e = copy.copy(self)
		e.transcript = str(self.transcript)
		return e

	def nextExon(self) :
		"""Returns the next exon of the transcript, or None if there is none"""
		try :
			return self.transcript.exons[self.number+1]
		except IndexError :
			return None
	
	def previousExon(self) :
		"""Returns the previous exon of the transcript, or None if there is none"""
		
		if self.number == 0 :
			return None
		
		try :
			return self.transcript.exons[self.number-1]
		except IndexError :
			return None
		
	def __str__(self) :
		return """EXON, id %s, number: %s, (start, end): (%s, %s), cds: (%s, %s) > %s""" %( self.id, self.number, self.start, self.end, self.CDS_start, self.CDS_end, str(self.transcript))

	def __len__(self) :
		return len(self.sequence)
