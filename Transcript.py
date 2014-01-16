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
	
	def __init__(self, importing = False) :
		
		if not importing :
			self.sequence = []
			self.CDNA = []
			self.UTR5 = []
			self.UTR3 = []
			prime5 = True
			for e in self.exons :
				self.sequence.append(e.sequence)
				if e.hasCDS() :
					self.CDNA.append(e.CDSSequence)
				else :
					if prime5 :
						self.UTR5.append(e.sequence)
					else :
						self.UTR3.append(e.sequence)					
			
			self.sequence = ''.join(self.sequence)
			self.CDNA = ''.join(self.CDNA)
			self.UTR5 = ''.join(self.UTR5)
			self.UTR3 = ''.join(self.UTR3)
			
			self.bin_CDNA = NucBinarySequence(self.CDNA)
			self.bin_sequence = NucBinarySequence(self.sequence)
			self.bin_UTR5 = NucBinarySequence(self.UTR5)
			self.bin_UTR3 = NucBinarySequence(self.UTR3)
		
			if len(self.CDNA) % 3 != 0 :
				self.flags = {'DUBIOUS' : True, 'CDNA_LEN_NOT_MULT_3': True} 
			else :
				self.flags = {'DUBIOUS' : False, 'CDNA_LEN_NOT_MULT_3': False} 
	
	def getNucleotideCodon(self, cdnaX1) :
		"Returns the entire codon of the nucleotide at pos cdnaX1 in the cdna, and the position of that nocleotide in the codon"
		return uf.getNucleotideCodon(self.CDNA, cdnaX1)

	def getCodon(self, i) :
		"returns the ith codon"
		return self.getNucleotideCodon(i*3)
	
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
		
	def __getitem__(self, i) :
		return self.sequence[i]
		
	def __len__(self) :
		return len(self.sequence)
	
	def __str__(self) :
		return """Transcript, id: %s, name: %s > %s""" %(self.id, self.name, str(self.gene))
		
