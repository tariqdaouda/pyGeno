from pyGenoObjectBases import *
from SNP import SNP_INDEL

import rabaDB.fields as rf
from tools import UsefulFunctions as uf
from tools.BinarySequence import NucBinarySequence

class Exon_Raba(pyGenoRabaObject) :
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

	_wrapped_class = Exon_Raba

	def __init__(self, *args, **kwargs) :
		"""An exon, the sequence is set according to gene strand, if it's '-' the sequence is the complement.
		A CDS is a couple of coordinates that lies inside of the exon.
		SNVsFilter is a fct that defines wich SNVs are included in the sequence.
		I expect [start, end[ (python) not something like [start, end](ensembl format)"""
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)
		self._load_sequencesTriggers = set(["UTR5", "UTR3", "CDS", "sequence"])

	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		if issubclass(objectType, SNP_INDEL) :
			f = RabaQuery(objectType, namespace = self._wrapped_class._raba_namespace)
			coolArgs['chromosomeNumber'] = self.chromosome.number
			coolArgs['start'] = self.start
			coolArgs['end'] = self.end
		
			if len(args) > 0 and type(args[0]) is types.ListType :
				for a in args[0] :
					if type(a) is types.DictType :
						f.addFilter(**a)
			else :
				f.addFilter(*args, **coolArgs)

			return f
		
		return pyGenoRabaObjectWrapper._makeLoadQuery(self, objectType, *args, **coolArgs)
	
	def _load_sequences(self) :
		seq = self.chromosome.sequence[self.start : self.end]
		if self.strand == '+' :
			self.sequence = seq
		else :
			self.sequence =  uf.reverseComplement(str(seq))
		
		if self.hasCDS() :
			start = self.CDS_start-self.start
			end = self.CDS_end-self.start
			
			if self.strand == '+' :
				self.UTR5 = self.sequence[:start]
				self.CDS = self.sequence[start:end]
				self.UTR3 = self.sequence[end:]
			else :
				self.UTR5 = self.sequence[:len(self.sequence)-end]
				self.CDS = self.sequence[len(self.sequence)-end:len(self.sequence)-start]
				self.UTR3 = self.sequence[len(self.sequence)-start:]
		else :
			self.UTR5 = ''
			self.CDS = ''
			self.UTR3 = ''

	def _load_bin_sequence(self) :
		self.bin_sequence = NucBinarySequence(self.sequence)
		self.bin_UTR5 =  NucBinarySequence(self.UTR5)
		self.bin_CDS =  NucBinarySequence(self.CDS)
		self.bin_UTR3 =  NucBinarySequence(self.UTR3)
		
	def hasCDS(self) :
		if self.CDS_start != None and self.CDS_end != None:
			return True
		return False

	def getCDSLength(self) :
		return len(self.CDS)

	def find(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.bin_sequence.find(sequence)

	def findAll(self, seqence):
		"""Returns a lits of all positions where sequence was found"""
		return self.bin_sequence.findAll(sequence)

	def findInCDS(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.bin_CDS.find(sequence)

	def findAllInCDS(self, seqence):
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
