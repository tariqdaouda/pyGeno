from tools import UsefulFunctions as uf
from tools.BinarySequence import AABinarySequence
import copy

class InvalidProtein(Exception):
	def __init__(self, protein, message):
		self.proteinId = protein.id
		self.protein = protein.sequence
		self.gene = protein.gene.symbol
		self.message = message
		self.binarySequence = AABinarySequence(self.sequence)
		
	def __str__(self):
		return """
		Description : %s
		protein_id : %s
		protein : %s
		gene_symbol : %s"""%(self.message, self.proteinId, self.protein, self.gene.symbol)

class Protein :
	def __init__(self, sequence, protId, transcript = None):

		self.sequence = sequence
		self.transcript = transcript
		self.id = protId
		self.binarySequence = AABinarySequence(self.sequence)
		self.updateBinarySequence = False
		
	def getSequence(self):
		return self.sequence
	
	def getDefaultSequence(self) :
		"""returns a version str sequence where only the last allele of each polymorphisms is shown"""
		return self.binarySequence.defaultSequence
	
	def getPolymorphisms(self) :
		return self.binarySequence.getPolymorphisms()

	def find(self, sequence):
		"""Returns the first occurence of sequence taking polymorphisme into account this is slower than the simple string search findString"""
		if self.updateBinarySequence :
			self.binarySequence = AABinarySequence(self.sequence)
			self.updateBinarySequence = False
		return self.binarySequence.find(sequence)

	def findAll(self, sequence):
		"""Returns all occurences of sequence taking polymorphisme into account this slower than the simple string search findStringAll"""
		if self.updateBinarySequence :
			self.binarySequence = AABinarySequence(self.sequence)
			self.updateBinarySequence = False
		return self.binarySequence.findAll(sequence)

	def findString(self, sequence) :
		"""return the first occurence of sequence using simple string search in sequence doesn't care about polymorphisme"""
		return self.sequence.find(sequence)

	def findStringAll(self, sequence):
		"""return all first occurences of sequence using simple string search in sequence doesn't care about polymorphisme"""
		return uf.findAll(self.sequence, sequence)
	
	def pluck(self) :
		"""Plucks the protein off the tree. Returns a protein identical to self but where the field .transcript has str(self.transcipt) as value,
		This makes the protein much more lighter in case you'd like to pickle it"""
		np = copy.copy(self)
		np.transcript  = str(self.transcript)
		return np
		
	def __getitem__(self, i) :
		return self.binarySequence.getChar(i)

	def __len__(self) :
		return len(self.binarySequence)

	def __str__(self) :
		return "Protein, id: %s -|- %s" %(self.id, str(self.transcript))
