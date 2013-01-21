import configuration as conf
from tools import UsefulFunctions as uf

from Chromosome import Chromosome
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
		
		self.genomes = {}
		for p in paths :
			self.genomes[p] = Genome(p, verbose)
			if not self.genomes[p].isLight :
				raise ValueError("Only light genomes can be mixed, %s is not light" % p)
				
	def empty(self) :
		self.genomes = {}
	
	def loadChromosome(self, number, dbSNPVersion = None, verbose = False) :
		for g in self.genomes.values():
			g.loadChromosome(number, dbSNPVersion, verbose)
	
	def loadAllChromosomes(self, dbSNPVersion, verbose = False) :
		for g in self.genomes.values():
			g.loadAllChromosomes(number, dbSNPVersion, verbose)

	def unloadChromosome(self, number) :
		for g in self.genomes.values():
			g.unloadChromosome(number)
	
	def getSequence(self, SNVsFilter, chroNumber, x1, x2 = None) :
		"""SNVsFilter is the function that defines the way the snv from different genomes must be mixed."""
		
		if x1 != None :
			if x2 == None :
				start, end = x1, x1 + 1
			elif x1 > x2 :
				start, end = x2, x1 
			else :
				start, end = x1, x2
				
			snps = {}
			for p in self.genomes :
				snps[p] = self.genomes.loadChromosome(chroNumber).getSNVsInRange(x1, x2)
				
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