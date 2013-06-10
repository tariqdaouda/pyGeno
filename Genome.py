import configuration as conf
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

class Genome :
	
	def __init__(self, path, reference = None, verbose = False) :
		"""path is a string that must have the following form: specie/genomeName ex: 'human/GRCh37.p2' or 'human/patient1'
		To know more about how to import new genomes please have a look at Importation.py. Beside if the genome name is
		'reference', ex : 'human/reference' the default reference genome will be loaded
		
		In case the genome is not complete (missing data for some chromosomes) or if your genome only contains a list of snps,
		(ex: sequencing results). To fill in the missing data pyGeno needs a complete genome to use as reference. As the specie
		has already been specified in path, you only need to pass the name of the reference genome ex: "GRCh37.p2".
		If no reference genome is specified and one is needed pyGeno will try to load the default reference genome for the specie
		as defined in pyGeno_SETTINGS['REFERENCE_GENOMES'][specie] (or more infos on defaults, see update_REFERENCE_GENOME() in configuration.py). 
		"""
		
		self.reset(path, reference, verbose)
	
	def reset(self, path, reference = None, verbose = False) :
		if verbose :
			print "Creating genome: " + path + "..."
		
		self.path = path
		self.specie = path.split('/')[0]
		self.name = path.split('/')[1]
		
		if self.name == 'reference' :
			self.name = conf.get_REFERENCE_GENOME(self.specie)
			
		if reference != None :
			self.reference = reference
		else :
			self.reference = conf.get_REFERENCE_GENOME(self.specie)
		
		self.absolutePath = conf.pyGeno_SETTINGS['DATA_PATH']+'/%s/genomes/%s' % (self.specie, self.name)
		self.referenceAbsolutePath = conf.pyGeno_SETTINGS['DATA_PATH']+'/%s/genomes/%s' % (self.specie, self.reference)	
		
		if os.path.exists(self.absolutePath + '/genomeChrPos.index') :
			f = open(self.absolutePath + '/genomeChrPos.index')
		elif os.path.exists(self.referenceAbsolutePath + '/genomeChrPos.index') :
			f = open(self.referenceAbsolutePath + '/genomeChrPos.index')
		else :
			raise GenomeError("Can't find genomeChrPos.index in neither %s nor the reference genome %s" % (path, self.reference), self.path)
		
		self.chrsData = {}
		startPosition = 0
		for l in f.readlines()[1:] :
			sl = l.split(';')
			self.chrsData[sl[0]] = ChrData_Struct(self, sl[0], sl[1], sl[2], sl[3])
			self.length = self.chrsData[sl[0]].x2
		
		f.close()
		
		self.empty()
	
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
		
	def getdbSNPPath(self) :
		return conf.pyGeno_SETTINGS['DATA_PATH']+'/%s/dbSNP/' % (self.specie)

	def empty(self) :
		self.chromosomes = {}
	
	def loadChromosome(self, numberStr, dbSNPVersion = None, verbose = False) :
		number = numberStr.upper()
		if number not in self.chromosomes.keys():
			if number != '' :
				self.chromosomes[number] = Chromosome(number, self, self.chrsData[number].x1, self.chrsData[number].x2, dbSNPVersion, verbose)

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
				print 'aaaa', c, chros
				return self.loadChromosome(c.number)
				
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
		return "Genome: %s, ref: %s" %(self.path, self.reference)
