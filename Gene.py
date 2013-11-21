import re, copy, sys, random

import configuration as conf

from rabaDB.setup import *
RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
from rabaDB.Raba import *
import rabaDB.fields as rf

from tools import UsefulFunctions as uf
from Protein import Protein
from Exon import Exon
from Transcript import Transcript

#from tools.BinarySequence import NucBinarySequence

class Gene(Raba) :
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	
	id = rf.PrimitiveField()
	name = rf.PrimitiveField()
	strand = rf.PrimitiveField()
	biotype = rf.PrimitiveField()
	
	genome = rf.RabaObjectField('Genome')
	chromosome = rf.RabaObjectField('Chromosome')
	transcripts = rf.RabaRelationField('Transcript')
	
	_raba_uniques = [('genome', 'id')]
	
	def __init__(self, *args, **fieldsSet) :
		Raba.__init__(self, **fieldsSet) 
	
	def loadTranscript(self, transId) :
		return self.transcripts[transId]
	
	def loadRandomTranscript(self) :
		"""returns a transcript at random"""
		k = int(random.random()*len(self.transcripts.keys()))
		key = self.transcripts.keys()[k]
		return self.loadTranscript(key)
	
	def getExons(self):
		return self.exons
		
	def loadProtein(self, protId) :
		"""Browse the transcript list looking for the one that code for the protein and returns
		a protein if found"""
		for t in self.transcripts.values() :
			if t.proteinId == protId :
				return t.loadProtein()
		return None 

	def getTranscriptIds(self) :
		return self.transcripts.keys()

	def getTranscripts(self, unDubiousOnly = False) :
		if not unDubiousOnly :
			return self.transcripts.values()
		return self.unDubiousTranscripts
		
	def getPolymorphisms(self) :
		ret = []
		for e in self.exons :
			ret.extend(e.getPolymorphisms())
		return ret

	def pluck(self) :
		"""Returns a plucked object. Plucks the gene off the tree, set the value of self.chromosome into str(self.chromosome). This effectively disconnects the object and
		makes it much more lighter in case you'd like to pickle it"""
		e = copy.copy(self)
		e.chromosome = str(self.chromosome)
		return e
	
	def __str__(self) :
		return "Gene, name: %s, id: %s, strand: %s > %s" %(self.name, self.id, self.strand, str(self.chromosome))
