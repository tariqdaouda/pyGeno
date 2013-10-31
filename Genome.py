import configuration as conf

from rabaDB.setup import *
RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
from rabaDB.Raba import *
import rabaDB.fields as rf

from tools import UsefulFunctions as uf

from Chromosome import Chromosome
#import sys, pickle, random, shutil, os, glob
import os, pickle
from tools import SingletonManager
from exceptions import *

class ChrData_Struct :
	
	def __init__(self, genome, number, x1, x2, length) :
		self. x1 = int(x1)
		self.x2 = int(x2)
		self.number = number
		self.length = int(length)
		
		#indexFp = conf.DATA_PATH+'%s/gene_sets/chr%s_gene_symbols.index.pickle'%(genome.getSpecie(), number)
		indexFp = conf.pyGeno_SETTINGS['DATA_PATH']+'%s/gene_sets/chr%s_gene_symbols.index.pickle'%(genome.specie, number)
		f = open(indexFp)
		if not SingletonManager.contains(indexFp) :
			SingletonManager.add(pickle.load(f), indexFp)	
		
		self.geneSymbolIndex = SingletonManager.get(indexFp)
		f.close()
	
	def hasGene(self, symbol) :
		return symbol in self.geneSymbolIndex.keys()

class Genome(Raba) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	
	id = rf.PrimitiveField()
	name = rf.PrimitiveField()
	specie = rf.PrimitiveField()
	reference = rf.PrimitiveField()
	chromosomes = rf.RabaListField()
	genomeSource = rf.PrimitiveField()
	packageInfos = rf.PrimitiveField()
	
	
	def __init__(self, *args, **fieldsSet) :
		Raba.__init__(self, **fieldsSet)
		
		self.verbose = False
		#self.absolutePath = conf.DATA_PATH+'/%s/genomes/%s' % (self.specie, self.name)
		#self.referencePath = conf.DATA_PATH+'/%s/genomes/reference' % (self.specie)
		
		#try :
		#	f = open(self.absolutePath + '/genomeChrPos.index')
		#except IOError:
		#	f = open(self.referencePath + '/genomeChrPos.index')
	
	def getChromosomesNumberList(self):
		return self.chrsData.keys()
	
	def getSequencePath(self) :
		return conf.pyGeno_SETTINGS['DATA_PATH']+'/%s/genomes/%s' % (self.specie, self.name)
	
	def getReferenceSequencePath(self) :
		return conf.pyGeno_SETTINGS['DATA_PATH']+'/%s/genomes/%s' % (self.specie, self.reference)
		
	def getGeneSetsPath(self) :
		return conf.pyGeno_SETTINGS['DATA_PATH']+'/%s/genomes/%s/gene_sets' % (self.specie, self.name)
	
	def getReferenceGeneSetsPath(self) :
		return conf.pyGeno_SETTINGS['DATA_PATH']+'/%s/genomes/%s/gene_sets' % (self.specie, self.reference)
		
	def loadChromosome(self, number, dbSNPVersion = None) :
		number = number.upper()
		if number not in self.chromosomes.keys():
			if number != '' :
				self.chromosomes[number] = Chromosome(number = number, genome = self, dbSNPVersion = dbSNPVersion)

		return self.chromosomes[number]
	
	def getChromosomes(self) :
		self.loadAllChromosomes()
		return self.chromosomes.values()
	
	def loadAllChromosomes(self, dbSNPVersion = None, verbose = False) :
		for c in self.chrsData.keys() :
			self.loadChromosome(c, dbSNPVersion, verbose)
	
	def unloadChromosome(self, number) :
		del(self.chromosomes[number])
	
	def getSequence(self, x1, x2) :
		c = self.chrsData[self.whereIsPosition(x1)]
		chro = self.loadChromosome(c.number)
		return chro.getSequence(x1 - c.x1, x2 - c.x1)
	
	def loadGene(self, symbol) :
		chro = whereIsGene(symbol)
		return chro.loadGene(symbol)
		
	def whereIsPosition(self, x1) :
		"""returns the nimber of the chormosome that correspond to pos x1"""
		for c in self.chrsData.values() :
			if c.x1 <= x1 <= c.x2 :
				return c.number
	
	def whereIsGene(self, geneSymbol, chromosomes = None) :
		""""Special request from Diana:
		Returns the chromosome that contains a given gene
		@param: geneSymbol a string
		@param: chromosomes list of chromosomes numbers
		if chromosomes is None, the whole list of genome's chromosome will be checked
		"""
		if chromosomes == None : 
			chros = self.chrsData.keys()
		else :
			chros = chromosomes 
		
		res = {}
		for c in chros :
			if self.chrsData[c].hasGene(geneSymbol) :
				return self.loadChromosome(c)
				
	def loadRandomChromosome(self, loadSnps = True) :
		"""Picks a random position in the genome and returns it's chromosome"""
		x1 = int(random.random()*self.length)
		return self.loadChromosome(self.whereIsPosition(x1), loadSnps)
		
	def isLoaded(self, number) :
		return number in self.chromosomes.keys()
		
	def __getitem__(self, chromo) :
		return self.chromosomes[chromo]

	def __len__(self) :
		"""Size of the genome in pb"""
		return self.length

	def __str__(self) :
		return "Genome: %s/%s" %(self.specie, self.name)
		#return "Genome: %s, ref: %s" %(self.path, self.reference)
