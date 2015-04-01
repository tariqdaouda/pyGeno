import configuration as conf

from pyGenoObjectBases import *
from SNP import SNP_INDEL

import rabaDB.fields as rf

from tools import UsefulFunctions as uf
from tools.BinarySequence import AABinarySequence
import copy

class Protein_Raba(pyGenoRabaObject) :
	"""The wrapped Raba object that really holds the data"""
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	id = rf.Primitive()
	name = rf.Primitive()

	genome = rf.RabaObject('Genome_Raba')
	chromosome = rf.RabaObject('Chromosome_Raba')
	gene = rf.RabaObject('Gene_Raba')
	transcript = rf.RabaObject('Transcript_Raba')

	def _curate(self) :
		if self.name != None :
			self.name = self.name.upper()

class Protein(pyGenoRabaObjectWrapper) :
	"""The wrapper for playing with Proteins"""
	
	_wrapped_class = Protein_Raba

	def __init__(self, *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)
		self._load_sequencesTriggers = set(["sequence"])

	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		if issubclass(objectType, SNP_INDEL) :
			f = RabaQuery(objectType, namespace = self._wrapped_class._raba_namespace)
			coolArgs['species'] = self.genome.species
			coolArgs['chromosomeNumber'] = self.chromosome.number
			coolArgs['start >='] = self.transcript.start
			coolArgs['start <'] = self.transcript.end
		
			if len(args) > 0 and type(args[0]) is types.ListType :
				for a in args[0] :
					if type(a) is types.DictType :
						f.addFilter(**a)
			else :
				f.addFilter(*args, **coolArgs)

			return f
		
		return pyGenoRabaObjectWrapper._makeLoadQuery(self, objectType, *args, **coolArgs)
	
	def _load_sequences(self) :
		self.sequence = uf.translateDNA(self.transcript.cDNA[:-3])
	
	def getSequence(self):
		return self.sequence

	def _load_bin_sequence(self) :
		self.bin_sequence = AABinarySequence(self.sequence)

	def getDefaultSequence(self) :
		"""Returns a version str sequence where only the last allele of each polymorphisms is shown"""
		return self.bin_sequence.defaultSequence

	def getPolymorphisms(self) :
		"""Returns a list of all polymorphisms contained in the protein"""
		return self.bin_sequence.getPolymorphisms()

	def find(self, sequence):
		"""Returns the position of the first occurence of sequence taking polymorphisms into account"""
		return self.bin_sequence.find(sequence)

	def findAll(self, sequence):
		"""Returns all the position of the occurences of sequence taking polymorphisms into accoun"""
		return self.bin_sequence.findAll(sequence)

	def findString(self, sequence) :
		"""Returns the first occurence of sequence using simple string search in sequence that doesn't care about polymorphisms"""
		return self.sequence.find(sequence)

	def findStringAll(self, sequence):
		"""Returns all first occurences of sequence using simple string search in sequence that doesn't care about polymorphisms"""
		return uf.findAll(self.sequence, sequence)

	def __getitem__(self, i) :
		return self.bin_sequence.getChar(i)
		
	def __len__(self) :
		return len(self.bin_sequence)

	def __str__(self) :
		return "Protein, id: %s > %s" %(self.id, str(self.transcript))
