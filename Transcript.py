import numpy as N
import re, copy, sys, random

import configuration as conf

from rabaDB.setup import *
RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
from rabaDB.Raba import *
import rabaDB.fields as rf

from tools import UsefulFunctions as uf
from Protein import Protein
from Exon import Exon

from tools.BinarySequence import NucBinarySequence


class Transcript(Raba) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	
	id = rf.Primitive()
	name = rf.Primitive()
	length = rf.Primitive()
	
	genome = rf.RabaObject('Genome')
	chromosome = rf.RabaObject('Chromosome')
	gene = rf.RabaObject('Gene')
	protein = rf.RabaObject('Protein')
	exons = rf.Relation('Exon')
	
	_raba_uniques = [('genome', 'id')]
	
	def __init__(self, *args, **fieldsSet) :
		pass
		"""self.gene = gene
		self.id = id
		self.name = name
		self.proteinId = proteinId
	
		self.exons = []
		self.codingExons = []
		self.sequence = ''
		self.CDNA = ''
		
		self.binCDNA = None
		self.binSequence = None
		
		self.firstCDSExon = None
		self.protein = None
		self.strand = self.gene.strand
		self.startCodon = None
		self.stopCodon = None
		
		#This dict contains a set au flags about the transcript
		#and his updated automaticly, and in theory should be modified directly.
		#the flag DUBIOUS is on if : LEN_NOT_MULT_3
		self.flags = {'DUBIOUS' : False, 'CDNA_LEN_NOT_MULT_3': False} 
		self.unclose()"""
	
	def appendExon(self, exon) :
		self.exons.append(exon)
		self.sequence += exon.sequence
		if exon.hasCDS() :
			if self.firstCDSExon == None :
				self.firstCDSExon = len(self.exons) -1
			seq = exon.getCDSSequence()
			if seq != '' :
				self.codingExons.append(exon)
				self.CDNA += seq

		if self.closed :
			self.__updateBinarySequences()
	
	"""def getFirstCdsNumber(self) :
		return self.exons[0].number 
	
	def getNumberOfExons(self) :
		return len(self.exons)
	
	def isProteinCoding(self):
		return (self.protein != '')
	"""
	
	def getRangeExons(self, dnaX1, dnaX2 = None) :
		"Returns the exons where dnaX1, dnaX2 is included. dnaX2 = None only dnaX1 is taken into account"
		ret = []
		if dnaX2 == None :
			for e in self.exons :
				if e.x1 <= dnaX1 and dnaX1 < e.x2 :
					ret.append(e)
		else:
			for e in self.exons :
				if (e.x1 <= dnaX1 and dnaX1 < e.x2) and (e.x1 <= dnaX2 and dnaX2 < e.x2) :
					ret.append(e)
		return ret
	
	def getNucleotideCodon(self, cdnaX1) :
		"Returns the entire codon of the nucleotide at pos cdnaX1 in the cdna, and the position of that nocleotide in the codon"
		return uf.getNucleotideCodon(self.CDNA, cdnaX1)

	def getCodon(self, codonNumber) :
		return self.getNucleotideCodon(codonNumber*3)
		
	def close(self) :
		"""After a transcript has been closed youcan still add exons but each new append would cause it
		recompute it's binary sequences, resulting in a overhead. A transcript that has not been closed
		cannot be searched using find and findAll. close() is autmaticly called by gene so shouldn't have to worry
		about it, but just in case you can unclose a transcript using unclose()"""
		self.__updateBinarySequences()
		self.closed = True
		self.__setFlags()
	
	def unclose(self) :
		self. __resetFlags()
		self.closed = False
	
	def __setFlags(self) :
		if len(self.CDNA)%3 != 0 :
			self.flags['CDNA_LEN_NOT_MULT_3'] = True
			self.flags['DUBIOUS'] = True
	
	def __resetFlags(self) :
		self.flags['CDNA_LEN_NOT_MULT_3'] = False
		self.flags['DUBIOUS'] = False
			
	def __updateBinarySequences(self) :
		if len(self.CDNA) > 0:
			self.binCDNA = NucBinarySequence(self.CDNA)
		self.binSequence = NucBinarySequence(self.sequence)
	
	def loadProtein(self) :
		if self.protein == '' :
			return None
		else :
			if self.protein == None :
				self.protein =  Protein(uf.translateDNA(self.CDNA), self.proteinId, self)
			return self.protein
	
	def find(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.binSequence.find(sequence)
	
	def findAll(self, seqence):
		"""Returns a lits of all positions where sequence was found"""
		return self.binSequence.findAll(sequence)

	def findInCDNA(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.binCDNA.find(sequence)
	
	def findAllInCDNA(self, seqence):
		"""Returns a lits of all positions where sequence was found"""
		return self.binCDNA.findAll(sequence)

	def getCDNALength(self):
		return len(self.CDNA)
		
	def CDNALenMultOf3(self) :
		return len(self.CDNA)%3 == 0
		
	def getCDNAStartPosition_todo(self):
		#TODO
		pass
		
	
	def getStartPosition_todo(self):
		#TODO
		pass
	
	def get5PrimeUtr_todo(self) :
		#TODO
		pass
		
	def get3PrimeUtr_todo(self) :
		#TODO
		pass
		
	def getFivePrimeUtrLength_todo(self):
		#TODO
		pass
		
	def getThreePrimeUtrLength_todo(self):
		#TODO
		pass
		
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
	
	def getExons(self):
		return self.exons
		
	def pluck(self):
		"""Returns a plucked object. Plucks the transcript off the tree, set the value of self.gene into str(self.gene). This effectively disconnects the object and
		makes it much more lighter in case you'd like to pickle it"""
		e = copy.copy(self)
		e.gene = str(self.gene)
		return e
	
	def getNbCodons(self) :
		return len(self.CDNA)/3
		
	def __getitem__(self, i) :
		return self.sequence[i]
		
	def __len__(self) :
		return len(self.sequence)
	
	def __str__(self) :
		return """Transcript, id: %s, name: %s > %s""" %(self.id, self.name, str(self.gene))
		
