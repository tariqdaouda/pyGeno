import re, random, sys,pickle
from tools.SegmentTree import SegmentTree as SegmentTree
from tools.SecureMmap import SecureMmap as SecureMmap
from tools import UsefulFunctions as uf
from tools import SingletonManager

import configuration as conf
from Gene import Gene
from SNP import *

from expyutils.GTFTools import GTFFile

class GeneNotFound(Exception):
	def __init__(self, chromosome, geneSymbol, message = ''):
		self.message = message
		self.symbol = geneSymbol
		self.chromosome = chromosome
		
	def __str__(self):
		return """
		Description : %s
		gene_symbol : %s
		chromosome : %s\n"""%(self.message, self.symbol, self.chromosome)

class RequestError(Exception):
	def __init__(self, message = ''):
		self.message = message
		
	def __str__(self):
		return """Request Error: %s"""%(self.message)

def defaultSNVsFilter(casavaSnp) :
	"""The default rule to decide wether to take the most probable genotype or the
	reference, always returns true"""
	return True

def defaultDbSNPsFilter(dbSnp) :
	"""The default rule to decide wether to take the most probable genotype or the
	snp, always returns true"""
	return True
	
class Chromosome :
	"""A class that represents a Chromosome
	Attention: private region support en retard par rapport au public"""

	def __init__(self, number, genome, dbSNPVersion = None, verbose = False) :
		self.reset(number, genome, dbSNPVersion, verbose)

	def emptyRegions(self) :
		"""Erases regions """
		self.SNPs = {}
		self.SNPData = None
		self.SNPIndexTree = None
		
	def empty(self) :
		"""Erases genes and regions"""
		self.dbSNPs = None
		self.genes = {}
		self.emptyRegions()
		
	def reset(self, number, genome, dbSNPVersion = None, verbose = False) :
		
		if verbose :
			print "Loading chromosome %s..."%(number)
		
		self.number = str(number).upper()
		self.genome = genome
		
		try :
			self.casavaSNPs = SNPFile('%s/chr%s.casavasnps'%(self.genome.getSequencePath(), self.number), CasavaSNP)
			refSeq = '%s/chr%s.dat'%(self.genome.getReferenceSequencePath(), self.number)
			
			if not SingletonManager.contains(refSeq) :
				SingletonManager.add(SecureMmap(refSeq), refSeq)
				
			self.data = SingletonManager.get(refSeq)
				
			self.isLight = True
		except IOError:
			try :
				self.data = SecureMmap('%s/chr%s.dat'%(self.genome.getSequencePath(), self.number))
			except IOError:
				print 'Warning : couldn\'t find local version of chromosome %s, loading reference instead...' % self.number
				self.data = SecureMmap('%s/chr%s.dat'%(self.genome.getReferenceSequencePath(), self.number))
				
			self.isLight = False

		f = open('%s/chr%s.gtf'%(self.genome.getGeneSetsPath(), self.number), 'r')

		self.gtfLines = f.readlines()
		f.close()
		
		self.geneSymbolIndex = self.genome.chrsData[self.number].geneSymbolIndex
		
		self.empty()
		
		if dbSNPVersion != None and dbSNPVersion != False and dbSNPVersion != '':
			self.loadSNPs(dbSNPVersion, verbose)
		elif verbose :
			print 'Not loading SNPs because told to...'
		
	def loadSNPs(self, version, verbose = False) :
		try :
			if verbose :
				print '--Loading SNPs...'
			snpFile = '%s/%s/chr%s.pygeno-dbSNP'%(self.genome.getdbSNPPath(), version,self.number)
			if not SingletonManager.contains(snpFile) :
				SingletonManager.add(SNPFile(snpFile, dbSNP), snpFile)		
			self.dbSNPs = SingletonManager.get(snpFile)
			
		except IOError :
			sys.stderr.write('Unable to load snp data for chr %s of genome %s' %(self.number, self.genome.name))
			sys.stderr.write('\t->Unable to resolve path: %s/chr%s.pygeno-dbSNP'%(self.genome.getdbSNPPath(), self.number))
			

	def hasGene(self, symbol) :
		return symbol in self.geneSymbolIndex.keys()
	
	def loadGene(self, symbol, SNVsFilter = None, verbose = False) :
		"""SNVsFilter is a fct tha takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros).
		If SNVsFilter != None it will automaticly trigger a complete reloading of the gene, even if it has already been loaded"""
		
		if symbol not in self.genes.keys() or SNVsFilter != None :
			try :
				l1, l2 = self.geneSymbolIndex[symbol].split(';')
				geneData = ''.join(self.gtfLines[int(l1) : int(l2)])
			except KeyError:
				raise GeneNotFound(self.number, symbol, 'Impossible to load Gene %s, not found in index' % symbol)
		
			gtf = GTFFile()
			gtf.parseStr(geneData)
			self.genes[symbol] = Gene(self, gtf, SNVsFilter , verbose)
		
		return self.genes[symbol]

	def unloadGene(self, symbol) :
		del(self.genes[symbol])
	
	def loadAllGenes(self, SNVsFilter = None, verbose = False) :
		"""SNVsFilter is a fct tha takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		for symbol in self.geneSymbolIndex.keys() :
			self.loadGene(symbol, SNVsFilter, verbose)

	
	def getGenes(self, SNVsFilter = None) :
		self.loadAllGenes(SNVsFilter)
		return self.genes.values()
	
	def getNucleotide(self, x1, SNVsFilter = None) :
		"""SNVsFilter is a fct that takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		if not self.isLight :
			return self.data[x1]
		
		if (SNVsFilter != None) :
			fct = SNVsFilter
		else :
			fct = defaultSNVsFilter
			
		snp = self.casavaSNPs.findSnp(x1)
		if snp != None :
			return snp['max_gt']
		return self.data[x1]

	def getSequence(self, x1, x2 = None, SNVsFilter = None) :
		"""SNVsFilter is a fct that takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		
		if x1 != None :
			if x2 == None :
				start, end = x1, x1 + 1
			elif x1 > x2 :
				start, end = x2, x1 
			else :
				start, end = x1, x2
				
			if not self.isLight :
				return self.data[start:end]
			else :
				if (SNVsFilter != None) :
					fct = SNVsFilter
				else :
					fct = defaultSNVsFilter
				
				snps = self.casavaSNPs.findSnpsInRange(start, end)
				data = None
				if snps != None :
					for snp in snps:
						if fct(snp) :
							if data == None :
								data = list(self.data[start:end])
							pos = snp['pos'] - start#-1
							data[pos] = snp['max_gt']

				if data != None :
					return ''.join(data)
				else :
					return self.data[start:end]
	
	def getSequence_dbSNP(self, x1, x2 = None, SNPsFilter = None) :
		"""SNPsFilter is a fct that takes a dbSNP SNP as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulDbSNPsFilter is used."""
		
		if x1 != None :
			if x2 == None :
				start, end = x1, x1 + 1
			elif x1 > x2 :
				start, end = x2, x1 
			else :
				start, end = x1, x2
				
			if self.dbSNPs == None :
				return self.data[start:end]
			else :
				if (SNPsFilter != None) :
					fct = SNPsFilter
				else :
					fct = defaultDbSNPsFilter
				
				snps = self.dbSNPs.findSnpsInRange(start, end)
				print snps
				data = None
				if snps != None :
					for snp in snps:
						if fct(snp) :
							if data == None :
								data = list(self.data[start:end])
							pos = snp['pos'] - start#-1
							print uf.getPolymorphicNucleotide(snp['alleles'])
							data[pos] = uf.getPolymorphicNucleotide(snp['alleles'])

				if data != None :
					return ''.join(data)
				else :
					return self.data[start:end]
					
	def getPolymorphismsInRange(self, x1, x2) :
		return self.casavaSNPs.findSnpsInRange(x1, x2)
		
	def getRandomCDS(self) :
		cds = CDS()
		found = False
		while not found :
			found = cds.setFromGTFLine(random.sample(self.gtfLines, 1)[0])
		
		return cds
	
	def getSNPsInRange(self, x1, x2 = None, dbSNPsFilter = None) :
		"""dbSNPsFilter is a fct that takes a dbSNP as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		res = []

		if self.dbSNPs == None :
			raise RequestError("No dbSNP database loaded, recrete the chromosome or use loadSNPs()")
			#return res
		else :
			if (dbSNPsFilter != None) :
				fct = dbSNPsFilter
			else :
				fct = defaultDbSNPsFilter

			snps = self.dbSNPs.findSnpsInRange(x1, x2)
			
			if snps != None :
				for snp in snps:
					if fct(snp) :
						res.append(snp)
		return res
	
	def getSNVsInRange(self, x1, x2 = None, SNVsFilter = None) :
		"""SNVsFilter is a fct that takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		res = []
			
		if not self.isLight :
			raise RequestError("Genome is light, there's no information on separate SNVs")
			#return res
		else :
			if (SNVsFilter != None) :
				fct = SNVsFilter
			else :
				fct = defaultSNVsFilter

			snps = self.casavaSNPs.findSnpsInRange(x1, x2)
			
			if snps != None :
				for snp in snps:
					if fct(snp) :
						res.append(snp)
		return res

		
	def findSequencePy(self, sequence) :
		return self.data.find(sequence)

	def getGeneIndex(self) :
		return self.geneSymbolIndex
	
	def loadRandomGene(self, SNVsFilter = None) :
		"""Loads a gene at random. SNVsFilter is a fct tha takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		
		k = int(random.random()*len(self.geneSymbolIndex.keys()))
		key = self.geneSymbolIndex.keys()[k]
		return self.loadGene(key, SNVsFilter)
		
	def __getitem__(self, i) :
		return self.genes[i]
		
	def __len__(self):
		return len(self.data)

	def __str__(self) :
		return "Chromosome: number %s -|- %s" %(self.number, str(self.genome))
