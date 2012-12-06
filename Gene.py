import numpy as N
import re, copy, sys, random

from tools import UsefulFunctions as uf
from Protein import Protein

from tools.BinarySequence import NucBinarySequence

#from SNP import GeneSNP
import PositionConverter as PC

class InvalidGene(Exception):
	def __init__(self, gene, message):
		self.message = message
		self.symbol = gene.symbol
		self.chromosome = gene.chromosome.number
		
	def __str__(self):
		return """
		Description : %s
		gene_symbol : %s
		chromosome : %s\n"""%(self.message, self.symbol, self.chromosome)

class InvalidExonComponent(Exception):
	def __init__(self, msg):
		self.msg = msg
		
	def __str__(self):
		return self.msg

class Exon:

	def __init__(self, x1, x2, number, transcript, type, SNVsFilter = None, startCodon = -1, stopCodon = -1) :
		r"""An exon, the sequence is set according to gene strand, if it's '-' the sequence is the complement.
		A CDS is a couple of coordinates that lies inside of the exon.
		SNVsFilter is a fct tha takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		
		self.type = type
		self.transcript = transcript
		self.startCodonPos = startCodon
		self.stopCodonPos = stopCodon
		self.x1 = int(x1)
		self.x2 = int(x2)
		self.number = int(number)
		self.CDS = None
		#self.threePrimeUtr = None
		#self.fivePrimeUtr = None

		seq = self.transcript.gene.chromosome.getSequence(x1, x2+1, SNVsFilter)
		if self.transcript.gene.strand == '+' :
			self.sequence = seq
		else :
			#self.sequence = uf.complementarySequence(seq[::-1])
			self.sequence = uf.reverseComplement(seq)
		
		#self.sequence = seq
	
	def getPolymorphisms(self) :
		return self.transcript.gene.chromosome.getPolymorphismsInRange(self.x1, self.x2)
	
	def hasCDS(self) :
		if self.CDS != None :
			return True
		return False
	
	def setCDS(self, x1, x2):
		self.CDS = (int(x1), int(x2))
		#if self.x1 != x1 :
		#	self.fivePrimeUtr = (self.x1, int(x1)-1)
		#if self.x2 != x2 :
		#	self.threePrimeUtr = (int(x2), self.x2)
		
	def getCDSLength(self) :
		return len(self.getCDSSequence())

	#def getFivePrimeUtrLength(self) :
	#	return len(self.getFivePrimeUtrSequence())
	
	#def getThreePrimeUtrLength(self) :
	#	return len(self.getThreePrimeUtrSequence())	
	
	def getCDSSequence(self) :
		try :
			if self.transcript.gene.strand == '+' :
				return self.sequence[self.CDS[0]-self.x1:self.CDS[1]-self.x1+1]
			else :
				x1 = self.CDS[0]-self.x1
				x2 = self.CDS[1]-self.x1
				return self.sequence[len(self.sequence)-1-x2:len(self.sequence)-x1]
		except TypeError:
			return ''
			#raise InvalidExonComponent("Invalid CDS, maybe the exon hasn't got any")

	#def getFivePrimeUtrSequence(self) :
	#	try :
	#		if self.transcript.gene.strand == '+' :
	#			return self.sequence[self.fivePrimeUtr[0]-self.x1:self.fivePrimeUtr[1]-self.x1+1]
	#		else :
	#			x1 = self.fivePrimeUtr[0]-self.x1
	#			x2 = self.fivePrimeUtr[1]-self.x1
	#			return self.sequence[len(self.sequence)-1-x2:len(self.sequence)-x1]
	#	except TypeError:
	#		return ''
	#		#raise InvalidExonComponent("Invalid fivePrimeUtr, maybe the exon hasn't got any")
	
	#def getThreePrimeUtrSequence(self) :
	#	try :
	#		if self.transcript.gene.strand == '+' :
	#			return self.sequence[self.threePrimeUtr[0]-self.x1:self.threePrimeUtr[1]-self.x1+1]
	#		else :
	#			x1 = self.threePrimeUtr[0]-self.x1
	#			x2 = self.threePrimeUtr[1]-self.x1
	#			return self.sequence[len(self.sequence)-1-x2:len(self.sequence)-x1]
	#	except TypeError:
	#		return ''
	#		#raise InvalidExonComponent("Invalid threePrimeUtr, maybe the exon hasn't got any")
	
	def __getitem__(self, i) :
		return self.sequence[i]
	
	def __str__(self) :
		return """EXON number: %d, x1: %d, x2: %d, cds: %s, transcript: %s""" %(self.number, self.x1, self.x2, str(self.CDS),self.transcript.name)

	def __len__(self) :
		return len(self.sequence)

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
	
	def __getitem__(self, i) :
		return self.sequence[i]
		
	def __str__(self) :
		return """Transcript in Chr: %s, Gene: %s, id: %s, name: %s""" %(self.gene.chromosome.number, self.gene.symbol, self.id, self.name)

	def __len__(self) :
		return len(self.sequence)
	
	def __repr__(self):
		return self.sequence
		
class Gene :
	
	def __init__(self, chromosome, gtfFile, SNVsFilter = None, verbose = False) :
		"""SNVsFilter is a fct tha takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		
		if verbose :
			print "Loading gene %s..."%(symbol)
		self.chromosome = chromosome
		
		#self.gtfFile = "".join(re.findall('.+"%s".+\n'%(symbol), gtfFile))		
		#print self.gtfFile
		#headerPattern = re.compile('([xXyY]|[0-9]{1,2})\s([^\s]+)\s.+([\+\-]).+"(ENS[^\;]+)".+gene_name "([^\;]+)"')
		#geneHeader = headerPattern.match(self.gtfFile)
		
		#fieldValue = re.compile('.*([^\s]+)\s"(["\s]+)".*')
		self.gtfFile = gtfFile
		self.type = self.gtfFile.get(0, 'coding_type')
		self.strand = self.gtfFile.get(0, 'strand')
		self.id = self.gtfFile.get(0, 'gene_id')
		self.symbol = self.gtfFile.get(0, 'gene_name')
		
		self.regionIndex = self.chromosome.getRegionIndex(self.symbol)
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
			#print '--', self.gtfFile[i]
			
			if transId not in self.transcripts.keys() :
				self.transcripts[transId] = Transcript(self, transId, transName, protId)
			else :
				if self.transcripts[transId].proteinId == '' :
					self.transcripts[transId].proteinId = protId
					
			if self.gtfFile.get(i, 'region_type') == 'exon' :
				#if (x1, x2) in self.exons.keys() :
				#	exon = self.exons[(x1, x2)]
				#else :
				exon = Exon(x1, x2, exonNumber, self.transcripts[transId], codingType, SNVsFilter)
				#self.exons[(x1, x2)] = exon
				self.exons.append(exon)
				#self.indexRegion(x1, x2, repr(exon), exon)
			else :
				for e in self.exons :
					if e.transcript.id == transId and e.x1 <= x1 and x2 <= e.x2:
						if self.gtfFile.get(i, 'region_type') == 'CDS':
							e.setCDS(x1, x2)
							#self.indexRegion(x1, x2, 'CDS of: %s' % (repr(e)), None)
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
	
	#def indexRegion(self, x1, x2, name = '', referencedObject = None):
	#	return self.regionIndex.addRegion(x1, x2, name, referencedObject)
	
	def loadTranscript(self, transId) :
		return self.transcripts[transId]
	
	def loadRandomTranscript(self) :
		"""returns a transcript at random"""
		k = int(random.random()*len(self.transcripts.keys()))
		key = self.transcripts.keys()[k]
		return self.loadTranscript(key)
	
	#def getTranscript(self, transId) :
	#	"""Just an alias for loadTranscript for the sack of back compatibility"""
	#	return self.getTanscripts(transId)
	
	#def getSNPs(self) :
	#	self.SNPs = self.regionIndex.getSNPs()
	#	return self.SNPs
	
	def getStartPosition(self) :
		return self.regionIndex.getX1()
	
	def getStopPosition(self) :
		return self.regionIndex.getX2()
	
	def __len__(self) :
		"""Sum of coding regions lengths"""
		return self.regionIndex.stree.getEffectiveLength()

	#def getCDNA(self) :
	#	"""Todo"""
	#	pass
		
	def loadProtein(self, protId) :
		"""Browse the transcript list looking for the one that code for the protein and returns
		a protein if found"""
		for t in self.transcripts.values() :
			#print t
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
