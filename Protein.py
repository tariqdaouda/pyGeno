import configuration as conf

from pyGenoObjectBases import *
from SNP import SNP_INDEL

import rabaDB.fields as rf

from tools import UsefulFunctions as uf
from tools.BinarySequence import AABinarySequence
import copy

class Protein_Raba(pyGenoRabaObject) :
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

	_wrapped_class = Protein_Raba

	def __init__(self, *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)
		self._load_sequencesTriggers = set(["sequence"])

	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		if issubclass(objectType, SNP_INDEL) :
			f = RabaQuery(objectType, namespace = self._wrapped_class._raba_namespace)
			coolArgs['chromosomeNumber'] = self.chromosome.number
			coolArgs['start'] = self.transcript.start
			coolArgs['end'] = self.transcript.end
		
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
		"""returns a version str sequence where only the last allele of each polymorphisms is shown"""
		return self.bin_sequence.defaultSequence

	def getPolymorphisms(self) :
		return self.bin_sequence.getPolymorphisms()

	def find(self, sequence):
		"""Returns the first occurence of sequence taking polymorphisme into account this is slower than the simple string search findString"""
		return self.bin_sequence.find(sequence)

	def findAll(self, sequence):
		"""Returns all occurences of sequence taking polymorphisme into account this slower than the simple string search findStringAll"""
		return self.bin_sequence.findAll(sequence)

	def findString(self, sequence) :
		"""return the first occurence of sequence using simple string search in sequence doesn't care about polymorphisme"""
		return self.sequence.find(sequence)

	def findStringAll(self, sequence):
		"""return all first occurences of sequence using simple string search in sequence doesn't care about polymorphisme"""
		return uf.findAll(self.sequence, sequence)

	def pluck(self) :
		"""Returns a plucked object. Plucks the protein off the tree, set the value of self.transcript into str(self.transcript). This effectively disconnects the object and
		makes it much more lighter in case you'd like to pickle it"""
		e = copy.copy(self)
		e.transcript = str(self.transcript)
		return e

	def __getitem__(self, i) :
		return self.bin_sequence.getChar(i)
		
	def __len__(self) :
		return len(self.bin_sequence)

	def __str__(self) :
		return "Protein, id: %s > %s" %(self.id, str(self.transcript))
