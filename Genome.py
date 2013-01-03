import configuration as conf
from tools import UsefulFunctions as uf

from Chromosome import Chromosome
import sys, pickle, random, shutil, os, glob

	
class ChromosomeNotFound(Exception) :
	def __init__(self, msg) :
		self.msg = msg
	def __str__(self) :
		return self.msg + '\n'
		
class ChrData_Struct :
	
	def __init__(self, genome, number, x1, x2, length) :
		self. x1 = int(x1)
		self.x2 = int(x2)
		self.number = number
		self.length = int(length)
	
		#f = open(conf.DATA_PATH+'/ensembl/%s/chr%s_gene_symbols.index.pickle'%(genome.getSpecie(), number))
		f = open(conf.DATA_PATH+'%s/gene_sets/chr%s_gene_symbols.index.pickle'%(genome.getSpecie(), number))
		self.geneSymbolIndex = pickle.load(f)
		f.close()
	
	def hasGene(self, symbol) :
		return symbol in self.geneSymbolIndex.keys()

class Genome :
	
	def __init__(self, path = 'human/reference', verbose = False) :
		self.reset(path, verbose)
		
	def reset(self, path = 'human/reference', verbose = False) :
		if verbose :
			print "Creating genome: " + path + "..."
		
		self.path = path
		self.specie = path.split('/')[0]
		self.name = path.split('/')[1]
		
		self.absolutePath = conf.DATA_PATH+'/%s/genomes/%s' % (self.specie, self.name)
		self.referencePath = conf.DATA_PATH+'/%s/genomes/reference' % (self.specie)
		#: WARNING ONLY FOR HUMANS!!!
		#self.chromosomeList = map(uf.intToStr, range(1, 23))
		#self.chromosomeList.extend(['x', 'y'])
		
		try :
			f = open(self.absolutePath + '/genomeChrPos.index')
		except IOError:
			f = open(self.referencePath + '/genomeChrPos.index')
			
		self.chrsData = {}
		for l in f.readlines()[1:] :
			sl = l.split(';')
			#print sl
			self.chrsData[sl[0]] = ChrData_Struct(self, sl[0], sl[1], sl[2], sl[3]) 
			self.length = self.chrsData[sl[0]].x2
		f.close()
		
		self.empty()
	
	def getChromosomesList(self):
		return self.chrsData.keys()
	
	def getSequencePath(self) :
		#return conf.DATA_PATH+'/ncbi/%s/sequences/%s' % (self.specie, self.name)
		return conf.DATA_PATH+'/%s/genomes/%s' % (self.specie, self.name)
	
	def getReferenceSequencePath(self) :
		#return conf.DATA_PATH+'/ncbi/%s/sequences/reference' % (self.specie)
		return conf.DATA_PATH+'/%s/genomes/reference' % (self.specie)
		
	def getGeneSetsPath(self) :
		#return conf.DATA_PATH+'/ensembl/%s' % self.getSpecie()
		return conf.DATA_PATH+'/%s/gene_sets' % self.getSpecie()
	
	def getdbSNPPath(self) :
		#return conf.DATA_PATH+'/ncbi/%s/dbSNP/' % (self.specie)
		return conf.DATA_PATH+'/%s/dbSNP/' % (self.specie)
		
	def getSpecie(self):
		return self.specie
		
	def empty(self) :
		self.chromosomes = {}
	
	def loadChromosome(self, number, dbSNPVersion = None, verbose = False) :
		if number not in self.chromosomes.keys():
			#try :
			if number != '' :
				self.chromosomes[number] = Chromosome(number, self, dbSNPVersion, verbose)
			#except IOError:
			#	raise ChromosomeNotFound('Unable to load chromosome :' + number)
		return self.chromosomes[number]
	
	def getChromosomes(self) :
		return self.chromosomes.values()
	
	def loadAllChromosomes(self, verbose = False) :
		for c in self.chrsData.keys() :
			self.loadChromosome(c, verbose)
	
	def unloadChromosome(self, number) :
		del(self.chromosomes[number])
	
	def getSequence(self, x1, x2) :
		ret = ''
		c = self.chrsData[self.whereIs(x1)]
		loaded = False	
		if not self.isLoaded(c.number) :
			loaded = True
		chr = self.loadChromosome(c.number)
		
		ret = chr.getSequence(x1 - c.x1, min(x2 - c.x1, len(chr)-1))
		
		if loaded :
			self.unloadChromosome(c.number)
		
		return ret
	
	def whereIs(self, x1) :
		"""returns the nimber of the chormosome that correspond to pos x1"""
		for c in self.chrsData.values() :
			if c.x1 <= x1 <= c.x2 :
				return c.number
	
	def whereIsGene(self, geneSymbol, chromosomes = None) :
		""""Special request from Diana:
		Returns a list of chormosomes indicating where the gene has been found
		@param: geneSymbol a string
		@param: chromosomes list of chromosomes numbers
		if chromosomes is None, the whole list of genome's chromosome will be checked
		"""
		return self.whereAreGenes([geneSymbol], chromosomes)[geneSymbol]
		
	def whereAreGenes(self, geneSymbols, chromosomes = None) :
		""""Special request from Diana:
		Returns a dictionary indicating where the genes have been found
		{Gene : chromosome list}
		@param: geneSymbols list of gene symbols
		@param: chromosomes list of chromosomes to look into (strings)
		if chromosomes is None, the whole list of genome's chromosome will be checked
		"""
		if chromosomes == None : 
			chros = self.chrsData.keys()
		else :
			chros = chromosomes 
		
		res = {}
		#wasLoaded = False
		for g in geneSymbols :
			res[g] = []
			for c in chros :
				if self.chrsData[c].hasGene(g) :
					res[g].append(c)
				"""if self.isLoaded(c) :
					if self.chromosomes[c].hasGene(g) :
						res[g].append(c)
				else :
					f = open('ensembl/%s/chr%s_gene_symbols.index.pickle'%(self.name, c))
					self.geneSymbolIndex = pickle.load(f)
					f.close()
					try :
						self.geneSymbolIndex[g]
						res[g].append(c)
					except KeyError:
						pass"""
						
		return res
	
	def loadRandomChromosome(self, loadSnps = True) :
		"""Picks a random position in the genome and returns it's chromosome"""
		x1 = int(random.random()*self.length)
		return self.loadChromosome(self.whereIs(x1), loadSnps)
		
	def isLoaded(self, number) :
		return number in self.chromosomes.keys()
		
	def __getitem__(self, chromo) :
		return self.chromosomes[chromo]

	def __len__(self) :
		"""Size of the genome in pb"""
		return self.length

	def __str__(self) :
		return "Genome: %s" %(self.path)