import numpy as N
import re, copy, sys, random

from tools import UsefulFunctions as uf
from Protein import Protein
from Exon import Exon

from tools.BinarySequence import NucBinarySequence


class Transcript :
	def __init__(self, gene, id, name, proteinId = '') :
		self.gene = gene
		self.id = id
		self.name = name
		self.proteinId = proteinId
	
		self.exons = []
		self.codingExons = []
		self.sequence = ''
		self.CDNA = ''
		#self.threePrimeUtr = ''
		#self.fivePrimeUtr = ''
		
		self.binCDNA = None
		self.binSequence = None
		
		self.firstCDSExon = None
		self.protein = None
		self.strand = self.gene.strand
		self.startCodon = None
		self.stopCodon = None
		
		#<7iyed>
		#this is set to None until getCodonUsage is called once
		self.codonUsage = None
		#</7iyed>
		self.lowAffinityCodonsNb = -1
		self.highAffinityCodonsNb = -1
		
		#This dict contains a set au flags about the transcript
		#and his updated automaticly, and in theory should be modified directly.
		#the flag DUBIOUS is on if : LEN_NOT_MULT_3
		self.flags = {'DUBIOUS' : False, 'CDNA_LEN_NOT_MULT_3': False} 
		self.unclose()
	
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
			
			#seq = exon.getFivePrimeUtrSequence()
			#if seq != '' :
			#	self.fivePrimeUtr += seq
			
			#seq = exon.getThreePrimeUtrSequence()
			#if seq != '' :
			#	self.threePrimeUtr += seq
		
		if self.closed :
			self.__updateBinarySequences()
	
	"""def getFirstCdsNumber(self) :
		return self.exons[0].number 
	
	def getNumberOfExons(self) :
		return len(self.exons)
	
	def isProteinCoding(self):
		return (self.protein != '')
	"""
	
	def getExons(self, dnaX1, dnaX2 = None) :
		"Returns the exons that the range dnaX1, dnaX2 belongs to. dnaX2 = None only dnaX1 is taken into account"
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
	
	def getCodon(self, cdnaX1) :
		"Returns the entire codon of the nucleotide at pos cdnaX1 in the cdna, and the position of that nocleotide in the codon"
		return uf.getCodon(self.CDNA, cdnaX1)

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
			#print "================", self.proteinId
			if self.protein == None :
				self.protein =  Protein(uf.translateDNA(self.CDNA), self.proteinId, self)
			return self.protein
	
	def find(self, sequence) :
		"""return the position of the first occurance of sequence"""
		#return self.sequence.find(sequence)
		return self.binSequence.find(sequence)
	
	def findAll(self, seqence):
		"""Returns a lits of all positions where sequence was found"""
		return self.binSequence.findAll(sequence)
		#return uf.findAll(sequence, self.sequence)

	def findInCDNA(self, sequence) :
		"""return the position of the first occurance of sequence"""
		return self.binCDNA.find(sequence)
	
	def findAllInCDNA(self, seqence):
		"""Returns a lits of all positions where sequence was found"""
		return self.binCDNA.findAll(sequence)

	def getCDNAStartPosition(self):
		""""Returns the start position of the coding sequence of the transcript"""
		return self.exons[self.firstCDSExon].CDS[0]
	
	def getStartPosition(self):
		""""Returns the start position of the transcript in the genome"""
		return self.exons[0].x1

	def getCDNALength(self):
		return len(self.CDNA)

	def get5PrimeUtr(self) :
		try :
			return self.sequence[:self.startCodon.x1 - self.exons[0].x1]
		except :
			return ''
	def get3PrimeUtr(self) :
		try :
			return self.sequence[self.stopCodon.x1 - self.exons[0].x1:]
		except :
			return ''	
	#def getFivePrimeUtrLength(self):
	#	return len(self.fivePrimeUtr)
		
	#def getThreePrimeUtrLength(self):
	#	return len(self.threePrimeUtr)
	
	#<7iyed>
	def getCodonAffinityMap(self, chunkRatio = 0.05) :
		chunks = []
		if len(self.CDNA) < 3 :
			return None
			
		#print len(self.CDNA)
		for i in range(int(1/chunkRatio)) :
			chunks.append({'low_aff':0, 'high_aff':0})
		
		efflen = len(self.CDNA)/3
		for i in range(0, len(self.CDNA), 3) :
			if len(self.CDNA[i:i+3]) == 3 :
				effI = i/3.
				c = int(N.floor(effI/(efflen*chunkRatio)))
				#print i, effI/(efflen*chunkRatio), effI, c, len(chunks), efflen, len(self.CDNA)
				if self.CDNA[i:i+3] in uf.codonTable.keys():
					if self.CDNA[i:i+3] in uf.lowAffinityCodons :
						chunks[c]['low_aff'] += 1
					else :
						chunks[c]['high_aff'] += 1
				else :
					for codon in uf.polymorphicCondonCombinaisons(self.CDNA[i:i+3]) :
						if codon in uf.lowAffinityCodons :
							chunks[c]['low_aff'] += 1
						else  :
							chunks[c]['high_aff'] += 1
		
		for i in range(len(chunks)) :
			if chunks[i]['high_aff']+chunks[i]['low_aff'] == 0 :
				chunks[i] = chunks[i-1]
	
		return chunks
		
	def getCodonAffinityMap_bck(self, chunkRatio = 0.05) :
		chunks = []
		if len(self.CDNA) < 3 :
			return None
			
		#print len(self.CDNA)
		for i in range(int(1/chunkRatio)) :
			chunks.append({'low_aff':0, 'high_aff':0})
		
		efflen = len(self.CDNA)/3
		for i in range(0, len(self.CDNA), 3) :
			if len(self.CDNA[i:i+3]) == 3 :
				effI = i/3.
				c = int(effI/(efflen*chunkRatio))#int(i/float(len(self.CDNA))*1/chunkRatio)
				#print i, i/3, (i/3)/(efflen*chunkRatio), c
				if self.CDNA[i:i+3] in uf.codonTable.keys():
					if self.CDNA[i:i+3] in uf.lowAffinityCodons :
						chunks[c]['low_aff'] += 1
					else :
						chunks[c]['high_aff'] += 1
				else :
					for codon in uf.polymorphicCondonCombinaisons(self.CDNA[i:i+3]) :
						if codon in uf.lowAffinityCodons :
							chunks[c]['low_aff'] += 1
						else  :
							chunks[c]['high_aff'] += 1
			#else:
			#	print self.CDNA[i:i+3]
		return chunks
		
	def getCodonUsage(self) :
		if self.codonUsage != None :
			return self.codonUsage
		
		for k in uf.codonTable.keys() :
			self.codonUsage[k] = 0
			
		for i in range(0, len(self.sequence), 3) :
			if self.sequence[i:i+3] in uf.codonTable.keys() :
				self.codonUsage[self.sequence[i:i+3]] += 1
				if self.sequence[i:i+3] in uf.lowAffinityCodons :
					self.lowAffinityCodonsNb += 1
				else :
					self.highAffinityCodonsNb += 1
			else :
				for c in uf.polymorphicCondonCombinaisons(self.sequence[i:i+3]) :
					self.codonUsage[c] += 1
					if c in uf.lowAffinityCodons :
						self.lowAffinityCodonsNb += 1
					else :
						self.highAffinityCodonsNb += 1
	
		return self.codonUsage
	#</7iyed>
	
	def pluck(self):
		"""Plucks the transcript off the tree. Returns a protein identical to self but where the field .gene has str(self.gene) as value,
		This makes the transcript much more lighter in case you'd like to pickle it"""
		nt = copy.copy(self)
		nt.gene = str(self.gene)
		return nt
		
	def __getitem__(self, i) :
		return self.sequence[i]
		
	def __len__(self) :
		return len(self.sequence)
	
	def __str__(self) :
		return """Transcript, id: %s, name: %s / %s""" %(self.id, self.name, str(self.gene))
		
