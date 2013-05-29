import numpy as N
import re, copy, sys, random

from tools import UsefulFunctions as uf
from Protein import Protein
from Exon import Exon
from Transcript import Transcript

from tools.BinarySequence import NucBinarySequence

class Gene :
	
	def __init__(self, chromosome, gtfFile, SNVsFilter = None, verbose = False) :
		"""SNVsFilter is a fct tha takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		
		if verbose :
			print "Loading gene %s..."%(symbol)
		self.chromosome = chromosome
		
		self.gtfFile = gtfFile
		self.type = self.gtfFile.get(0, 'coding_type')
		self.strand = self.gtfFile.get(0, 'strand')
		self.id = self.gtfFile.get(0, 'gene_id')
		self.symbol = self.gtfFile.get(0, 'gene_name')
		
		self.SNPs = {}
		self.exons = []
		self.transcripts = {}
		self.unDubiousTranscripts = []
		
		for i in range(len(self.gtfFile)) :
			x1 = int(self.gtfFile.get(i, 'x1')) -1
			x2 = int(self.gtfFile.get(i, 'x2')) -1
			transId = self.gtfFile.get(i, 'transcript_id')
			transName = self.gtfFile.get(i, 'transcript_name')
			exonNumber = self.gtfFile.get(i, 'exon_number')
			protId = self.gtfFile.get(i, 'protein_id')
			codingType = self.gtfFile.get(i, 'coding_type')
			
			if transId not in self.transcripts.keys() :
				self.transcripts[transId] = Transcript(self, transId, transName, protId)
			else :
				if self.transcripts[transId].proteinId == '' :
					self.transcripts[transId].proteinId = protId
					
			if self.gtfFile.get(i, 'region_type') == 'exon' :
				exon = Exon(x1, x2, exonNumber, self.transcripts[transId], codingType, SNVsFilter)
				self.exons.append(exon)
			else :
				for e in self.exons :
					if e.transcript.id == transId and e.x1 <= x1 and x2 <= e.x2:
						if self.gtfFile.get(i, 'region_type') == 'CDS':
							e.setCDS(x1, x2)
						elif self.gtfFile.get(i, 'region_type') == 'start_codon':
							e.startCodon = x1
							e.transcript.startCodon = e
						elif self.gtfFile.get(i, 'region_type') == 'stop_codon':
							e.stopCodon = x1
							e.transcript.stopCodon = e
		if verbose :
			print "Appending exons to transcripts of gene : %s..."%(symbol)
		
		for e in self.exons:
			self.transcripts[e.transcript.id].appendExon(e)
	
		for tid in self.transcripts :
			if verbose :
				print "closing transcript : %s..."%(tid)
			self.transcripts[tid].close()
			if not self.transcripts[tid].flags['DUBIOUS'] :
				self.unDubiousTranscripts.append(self.transcripts[tid])
	
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
		return "Gene, symbol: %s, id: %s, strand: %s / %s" %(self.symbol, self.id, self.strand, str(self.chromosome))
