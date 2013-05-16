import os, sys, random
from types import *

import configuration as conf
from Gene import Gene
from SNP import *
from exceptions import *

from tools.SegmentTree import SegmentTree as SegmentTree
from tools.SecureMmap import SecureMmap as SecureMmap
from tools import UsefulFunctions as uf
from tools import SingletonManager
from tools.GTFTools import GTFFile


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

	def __init__(self, number, genome, x1, x2, dbSNPVersion = None, verbose = False) :
		"""x1, x2 are the prosition of the chromosome in the genome"""
		self.reset(number, genome, x1, x2, dbSNPVersion, verbose)
		
	def reset(self, number, genome, x1, x2, dbSNPVersion = None, verbose = False) :
		"""x1, x2 are the prosition of the chromosome in the genome"""
		if verbose :
			print "Loading chromosome %s..."%(number)
		
		self.number = str(number).upper()
		self.genome = genome
		self.x1 = int(x1)
		self.x2 = int(x2)
		self.genes = {}
		
		self.__loadSequence()
		
		gtfFp = self.__getAnnotations()
		f = open(gtfFp, 'r')
		if not SingletonManager.contains(gtfFp) :
			SingletonManager.add(f.readlines(), gtfFp)	
		self.gtfLines = SingletonManager.get(gtfFp)
		f.close()
		
		self.geneSymbolIndex = self.genome.chrsData[self.number].geneSymbolIndex

		if dbSNPVersion != None and dbSNPVersion != False and dbSNPVersion != '':
			self.loadSNPs(dbSNPVersion, verbose)
		elif verbose :
			print 'Not loading SNPs because told to...'
	
	def __loadSequence(self):
		if os.path.exists('%s/chr%s.dat'%(self.genome.getSequencePath(), self.number)) :
			self.data = self.__getHeavySequence('%s/chr%s.dat'%(self.genome.getSequencePath(), self.number))
			if self.data != None :
				self.isLight = False
				return True
		elif os.path.exists('%s/chr%s.casavasnps'%(self.genome.getSequencePath(), self.number)) :
			self.casavaSNPs = SNPFile('%s/chr%s.casavasnps'%(self.genome.getSequencePath(), self.number), CasavaSNP)
			self.isLight = True
		else :
			sys.stderr.write('Warning : couldn\'t find local version of chromosome %s, attempting load reference instead...' % self.number)
		
		self.data = self.__getHeavySequence('%s/chr%s.dat'%(self.genome.getReferenceSequencePath(), self.number))
		if self.data == None :
			raise ChromosomeError("Unable to load chromosome %s! Impossible to find reference sequence" % self.number, self.number)
		
		self.isLight = False
		return True
	
	def __getHeavySequence(self, path) :
		try :
			if not SingletonManager.contains(path) :
				SingletonManager.add(SecureMmap(path), path)
				
			return SingletonManager.get(path)
		except :
			return None
	
	def __getAnnotations(self) :
		p = '%s/chr%s.gtf'%(self.genome.getGeneSetsPath(), self.number)
		rp = '%s/chr%s.gtf'%(self.genome.getReferenceGeneSetsPath(), self.number)
		if os.path.exists(p) :
			return p
		elif os.path.exists(rp) :
			return rp
		
		raise ChromosomeError("Unable to load chromosome %s! Can't find gene annotion sequence neither in %s or %s" % (self.number, p, rp), self.number)
		
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
		
		assert type(x1) is IntType
		assert type(x2) is IntType
		
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
							snp['max_gt'] = uf.getPolymorphicNucleotide(snp['max_gt'])
							data[pos] = snp['max_gt']

				if data != None :
					return ''.join(data)
				else :
					return self.data[start:end]
	
	def getSequence_dbSNP(self, x1, x2 = None, SNPsFilter = None) :
		"""SNPsFilter is a fct that takes a dbSNP SNP as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulDbSNPsFilter is used."""
		
		assert type(x1) is IntType
		assert type(x2) is IntType
		
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
		assert type(x1) is IntType
		assert type(x2) is IntType
		
		res = []

		if self.dbSNPs == None :
			raise RequestError("No dbSNP database loaded through loadSNPs()")
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
		assert type(x1) is IntType
		assert type(x2) is IntType
		
		res = []
			
		if not self.isLight :
			raise RequestError("Genome is light, there's no information on separate SNVs")
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
	
	def pluck(self) :
		"""Returns a plucked object. Plucks the chromosome off the tree, set the value of self.genome into str(self.genome). This effectively disconnects the object and
		makes it much more lighter in case you'd like to pickle it"""
		e = copy.copy(self)
		e.genome = str(self.genome)
		return e
		
	def __getitem__(self, i) :
		return self.genes[i]
		
	def __len__(self):
		return len(self.data)

	def __str__(self) :
		return "Chromosome: number %s / %s" %(self.number, str(self.genome))
