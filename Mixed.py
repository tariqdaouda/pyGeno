import configuration as conf
from tools import UsefulFunctions as uf

from Chromosome import *
from Genome import *
import sys, pickle, random, shutil, os, glob

def defaultMixedSNVFilter(snvDict) :
	"""Default function returns a concatenation of all snvs"""
	snvs = {}
	for p in snvsDict:
		for snv in snvsDict[p] :
			if snv['Qmax_gt'] >= 20 :
				if snv['pos'] not in snvs :
					snvs[snv['pos']] = snv
				else :
					snvs[snv['pos']]['max_gt'] = uf.getPolymorphicNucleotide(snvs[snv['pos']]['max_gt']+snv['max_gt'])
	return snvs.values()
	
class MixedGenome :
	"""With this class you can mix serveral light genomes an tune the getSequence to set the way the polymorphims
	of different genomes are mixed"""
	
	def __init__(self, paths, verbose = False) :
		self.reset(paths, verbose)
		
	def reset(self, paths, verbose = False) :
		
		if len(paths) < 1 :
			raise ValueError("No genome path defined")
		
		if verbose :
			print "Loading mixed genomes"
	
		self.empty()
		self.specie = None
		for p in paths :
			if verbose :
				print 'loading Genome', p
				
			self.genomes[p] = Genome(p, verbose)
			if self.specie == None :
				self.specie = self.genomes[p].specie
			elif self.specie != self.genomes[p].specie :
				raise ValueError('All genomes must be of the same specie %s != %s (set verbose to true for more info)' % (self.specie, self.genomes[p].specie))
			
		self.chromosomesNumberList = self.genomes[paths[0]].getChromosomesNumberList()
	
	def empty(self) :
		self.genomes = {}
		self.chromosomes = {}
	
	def getGenomePaths(self):
		return self.genomes.keys()
		
	def loadChromosome(self, number, dbSNPVersion = None, verbose = False) :
		self.chromosomes[number] = MixedChromosome(number, self, dbSNPVersion, verbose)
		return self.chromosomes[number]
		
	def loadAllChromosomes(self, dbSNPVersion = None, verbose = False) :
		for number in self.chromosomesNumberList :
			self.chromosomes[number] = MixedChromosome(number, self, dbSNPVersion, verbose)

	def unloadChromosome(self, number) :
		for g in self.genomes.values():
			g.unloadChromosome(number)

	def getChromosomesNumberList(self) :
		return self.chromosomesNumberList
	
	def getChromosomes(self):
		return self.chromosomes.values()
		
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
		self.genes = {}
		
		if verbose :
			print "Loading mixed chromosomes", number
	
		for p in self.mixedGenomes.getGenomePaths() :	
			#if verbose :
			#	print '\t for genome', p
			self.chromosomes[p] = self.mixedGenomes[p].loadChromosome(number, dbSNPVersion, verbose)
			if not self.chromosomes[p].isLight :
				raise ValueError("Only light chromosomes can be mixed, chromosomes %s of genome %s is not light" % (number, p))
		
		tmpChro = self.chromosomes.values()[0]
		self.geneSymbolIndex = tmpChro.geneSymbolIndex
		self.length = len(tmpChro)
		self.x1 = tmpChro.x1
		self.x2 = tmpChro.x2
		self.gtfLines = tmpChro.gtfLines
		self.data = tmpChro.data
		
	def loadGene(self, symbol, SNVsFilter = None, verbose = False) :
		"""SNVsFilter is ftc that takes a dictionnary of SNVs : genomePath => list of snvs, and then returns a list of the
		selected snvs. If None, the function defaultMixedSNVFilter will be used"""
		
		if SNVsFilter == None :
			fct = defaultMixedSNVFilter
		else :
			fct = SNVsFilter
			
		#if symbol not in self.genes.keys() and SNVsFilter == None :
		#	raise ValueError('Gene %s has not been loaded before, please specify a SNVsFilter function' % symbol)
		
		if symbol not in self.genes.keys():
			try :
				l1, l2 = self.geneSymbolIndex[symbol].split(';')
				geneData = ''.join(self.gtfLines[int(l1) : int(l2)])
			except KeyError:
				raise GeneNotFound(self.number, symbol, 'Impossible to load Gene %s, not found in index' % symbol)
		
			gtf = GTFFile()
			gtf.parseStr(geneData)
			self.genes[symbol] = Gene(self, gtf, fct , verbose)
		
		return self.genes[symbol]
	
	def loadAllGenes(self, SNVsFilter = None, verbose = False ):
		"""SNVsFilter is ftc that takes a dictionnary of SNVs : genomePath => list of snvs, and then returns a list of the
		selected snvs. If None, the function defaultMixedSNVFilter will be used"""
		for symbol in self.geneSymbolIndex.keys() :
			self.loadGene(symbol, SNVsFilter, verbose)
	
	def getGenes(self) :
		return self.genes.values()
		
	def getSNVsInRange(self, x1, x2, SNVsFilter = None) :
		"""SNVsFilter is ftc that takes a dictionnary of SNVs : genomePath => list of snvs, and then returns a list of the
		selected snvs. If None, the function defaultMixedSNVFilter will be used"""
		
		if SNVsFilter == None :
			fct = defaultMixedSNVFilter
		else :
			fct = SNVsFilter
			
		if x2 == None :
			start, end = x1, x1 + 1
		elif x1 > x2 :
			start, end = x2, x1 
		else :
			start, end = x1, x2
		
		individualSNPS = {}
		for p in self.chromosomes:
			individualSNPS[p] = self.chromosomes[p].getSNVsInRange(start, end)
		
		return fct(individualSNPS)
		
	def getSequence(self, x1, x2, SNVsFilter = None) :
		"""SNVsFilter is ftc that takes a dictionnary of SNVs : genomePath => list of snvs, and then returns a list of the
		selected snvs. If None, the function defaultMixedSNVFilter will be used"""
		
		if x2 == None :
			start, end = x1, x1 + 1
		elif x1 > x2 :
			start, end = x2, x1 
		else :
			start, end = x1, x2
		
		if SNVsFilter == None :
			fct = defaultMixedSNVFilter
		else :
			fct = SNVsFilter
			
		finalSNPS = self.getSNVsInRange(start, end, fct)
		if len(finalSNPS) > 0 :
			data = list(self.data[start:end])
			for snp in finalSNPS:
				#if len(snp['max_gt']) != 1 :
				#	raise ValueError("len(snp['max_gt']) != 1 (= %s), it should be one single nucleotide (single or polymorphic)" % snp['max_gt'])
				snp['max_gt'] = uf.getPolymorphicNucleotide(snp['max_gt'])
				pos = snp['pos'] - start#-1
				data[pos] = snp['max_gt']
			return ''.join(data)
		
		return self.data[start:end]

	def __len__(self) :
		return self.length
		
	def __str__(self) :
		return "Mixed Chromosome: %s" %(self.number)
