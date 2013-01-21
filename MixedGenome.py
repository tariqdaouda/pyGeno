import configuration as conf
from tools import UsefulFunctions as uf

from Chromosome import *
from Genome import *
import sys, pickle, random, shutil, os, glob
	
class MixedGenome :
	"""With this class you can mix serveral light genomes an tune the getSequence to set the way the polymorphims
	of different genomes are mixed"""
	
	def __init__(self, paths, verbose = False) :
		self.reset(paths, verbose)
		
	def reset(self, paths, verbose = False) :
		
		if len(paths) < 1 :
			raise ValueError("No genome path defined")
		
		if verbose ;
			print "Loading mixed genomes"
	
		self.genomes = {}
		for p in paths :
			if verbose :
				print '\t Genome', p
			self.genomes[p] = Genome(p, verbose)
			if not self.genomes[p].isLight :
				raise ValueError("Only light genomes can be mixed, %s is not light" % p)
				
	def empty(self) :
		self.genomes = {}
		self.chromosomes = {}
	
	def loadChromosome(self, number, dbSNPVersion = None, verbose = False) :
		self.chromosomes[number] = MixedChromosome(number, self, dbSNPVersion, verbose)
		return self.chromosomes[number]
		
	def loadAllChromosomes(self, dbSNPVersion = None, verbose = False) :
		for number in gself.genomes.values().getChromosomesList
			self.chromosomes[number] = MixedChromosome(number, self, dbSNPVersion, verbose)

	def unloadChromosome(self, number) :
		for g in self.genomes.values():
			g.unloadChromosome(number)

	def __getitem__(self, val) :
		return self.genomes[val]

	def __len__(self) :
		"""Returns the number of loaded genomes"""
		return len(self.genomes)

	def __str__(self) :
		s = []
		for g in self.genomes.values():
			s.append('\t' + str(g))
		
		s = '\n'.join(s)
		return "Mixed Genomes:\n%s" %(s)

class MixedChromosome :
	
	def __init__(self, number, mixedGenomes, dbSNPVersion = None, verbose = False) :
		
		self.number = number
		self.mixedGenomes = mixedGenomes
		self.chromosomes = {}
	
		if verbose ;
			print "Loading mixed chromosomes", number
	
		for p in self.mixedGenomes :	
			if verbose :
				print '\t for genome', p
			self.chromosomes[p] = self.mixedGenomes[p].loadChromosome(number, self.mixedGenomes[p], dbSNPVersion, verbose)
		
		tmpChro = self.chromosomes.values()[0]
		self.geneSymbolIndex = tmpChro.geneSymbolIndex
		self.genes = {}
		
	def loadGene(self, symbol, SNVsFilter, verbose = False) :
		"""SNVsFilter defines how the SNVs of the different genomes are mixed together. It takes an dictionary of snvs 
		for each chromosome in the form: snvs[genomePath] => listOfSNVs (the list that you would get if you called getSNVsInRange on chromosome),
		an returns a list of the selected SNVs to be applied"""
		
		if symbol not in self.genes.keys():
			try :
				l1, l2 = self.geneSymbolIndex[symbol].split(';')
				geneData = ''.join(self.gtfLines[int(l1) : int(l2)])
			except KeyError:
				raise GeneNotFound(self.number, symbol, 'Impossible to load Gene %s, not found in index' % symbol)
		
			gtf = GTFFile()
			gtf.parseStr(geneData)
			self.genes[symbol] = Gene(self, gtf, SNVsFilter , verbose)
		
		return self.genes[symbol]
	
	def getSequence(self, SNVsFilter, x1, x2 = None) :
		"""SNVsFilter defines how the SNVs of the different genomes are mixed together. It takes an dictionary of snvs 
		for each chromosome in the form: snvs[genomePath] => listOfSNVs (the list that you would get if you called getSNVsInRange on chromosome),
		an returns a list of the selected SNVs to be applied"""
		
		if x2 == None :
			start, end = x1, x1 + 1
		elif x1 > x2 :
			start, end = x2, x1 
		else :
			start, end = x1, x2
		
		individualSNPS = {}
		for p in self.chromosomes:
			individualSNPS[p] = self.chromosomes.getSNPsInRange(start, end)
		
		finalSNPS = SNVsFilter(individualSNPS)
		if len(finalSNPS) > 0 :
			data = list(self.data[start:end])
			for snp in finalSNPS:
				pos = snp['pos'] - start#-1
				data[pos] = snp['max_gt']
			return ''.join(data)
		
		return self.data[start:end]
	
	def __str__(self) :
		return 'Mixed chromosomes number: %s\n %s' %(self.number, str(self.mixedGenomese))