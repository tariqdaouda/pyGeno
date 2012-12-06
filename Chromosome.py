import re, random, sys,pickle
from tools.SegmentTree import SegmentTree as SegmentTree
from tools.SecureMmap import SecureMmap as SecureMmap
from tools import UsefulFunctions as uf
from tools import SingletonManager as SM

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
		
#import xml.etree.cElementTree as ET
def defaultSNVsFilter(casavaSnp) :
	"""The default rule to decide wether to take the most probable genotype or the
	reference, always returns true"""
	return True

def defaultSNPsFilter(dbSnp) :
	"""The default rule to decide wether to take the most probable genotype or the
	snp, always returns true"""
	return True
	
class Chromosome :
	"""A class that represents a Chromosome
	Attention: private region support en retard par rapport au public"""

	def __init__(self, number, genome, loadSNPs = True, verbose = False) :
		self.reset(number, genome, loadSNPs, verbose)

	def emptyRegions(self) :
		"""Erases regions """
		self.SNPs = {}
		self.SNPData = None
		self.SNPIndexTree = None
		
	def empty(self) :
		"""Erases genes and regions"""
		self.genes = {}
		self.emptyRegions()
		
	def reset(self, number, genome, loadSNPs = True, verbose = False) :
		
		if verbose :
			print "Loading chromosome %s..."%(number)
		
		self.number = str(number).upper()
		self.genome = genome
	
		#self.number = self.number.replace('_alt', '')
		
		try :
			#self.casavaSNPs = CasavaSNPs('%s/chr%s.casavasnps'%(self.genome.getSequencePath(), self.number))
			self.casavaSNPs = SNPFile('%s/chr%s.casavasnps'%(self.genome.getSequencePath(), self.number), CasavaSNP)
			refSeq = '%s/chr%s.dat'%(self.genome.getReferenceSequencePath(), self.number)
			self.data = SM.get(refSeq)
			if self.data == None :
				SM.add(SecureMmap(refSeq), refSeq)
				self.data = SM.get(refSeq)
				
			self.isLight = True
			#print "aaa", '%s/chr%s.dat'%(self.genome.getReferenceSequencePath(), self.number), len(self.data )
		except IOError:
			try :
				self.data = SecureMmap('%s/chr%s.dat'%(self.genome.getSequencePath(), self.number))
			except IOError:
				print 'Warning : couldn\'t find local version of chromosome %s, loading reference instead...' % self.number
				self.data = SecureMmap('%s/chr%s.dat'%(self.genome.getReferenceSequencePath(), self.number))
				
			self.isLight = False
			
		#f = open('pyGenoData/ncbi/%s/sequences//chr%s.dat'%(self.genome.name, number), 'r+b')
		#f = open('%s/chr%s.dat'%(self.genome.getSequencePath(), number), 'r+b')
		#mmap.mmap(f.fileno(), 0)#f.read()
		#f.close()
		
		#f = open('pyGenoData/ensembl/%s/chr%s.gtf'%(self.genome.name, number), 'r+b')
		#f = open('%s/chr%s.gtf'%(self.genome.getAnnotationPath(), self.number), 'r')
		f = open('%s/chr%s.gtf'%(self.genome.getGeneSetsPath(), self.number), 'r')
		#self.gtf = mmap.mmap(f.fileno(), 0)
		self.gtfLines = f.readlines()
		f.close()
		
		#f = open('ensembl/%s/chr%s_gene_symbols.index.pickle'%(self.genome.name, number))
		#self.geneSymbolIndex = pickle.load(f)
		#f.close()
		self.geneSymbolIndex = self.genome.chrsData[self.number].geneSymbolIndex
		
		self.empty()
		if loadSNPs :
			self.loadSNPs(verbose)
		elif verbose :
			print 'Not loading SNPs because told to...'
		
	def loadSNPs(self, verbose = False) :
		if self.SNPData == None :
			try :
				if verbose :
					print '--Loading SNPs...'
				
				specie = self.genome.getSpecie()
				
				f = open('%s/chr%s.sort.pygeno-dbSNP'%(self.genome.getSNPsPath(), self.number))
				self.SNPData = f.readlines()
				f.close()
				if verbose :
					print '--Loading SNP Index...'
				self.SNPIndexTree = SNPIndexTree('%s/chr%s.sort.pygenosnp.index' %(self.genome.getSNPsPath(), self.number))
			except IOError:
				print self.genome.getSNPsPath()
				sys.stderr.write('Unable to load snp data for chr %s of genome %s' %(self.number, self.genome.name))

	#def materialiseSNP(self, line) :
	#	"""creates a SNP object from data line. The snp is memorised
	#	by the chromosome
	#	@return the SNP object"""
	#	#print self.SNPData[line]
	#	sl = self.SNPData[line].split('#v#')
	#	snpL = sl[0].split(';')
	#	#'chr;chr-pos;rsid;alleles;orient;assembly;validated'
	#	snp = SNP(self.number, snpL[1], snpL[2], snpL[3], snpL[4], snpL[5], snpL[6])
	#	for lv in sl[1:] :
	#		slv = lv.split(';')
	#		#print 'a', slv
	#		#'gene-symbol;fxn-class;allele;residue;aa_position;frame'
	#		#(snp, geneSymbol, function, allele, residue, aaPosition, frame)
	#		snp.addVariant(slv[0], slv[1], slv[2], slv[3], slv[4], slv[5])
	#		
	#	self.SNPs[snpL[1]] = snp
	#	
	#	return snp
	
	def hasGene(self, symbol) :
		return symbol in self.geneSymbolIndex.keys()
	
	def loadGene(self, symbol, SNVsFilter = None, verbose = False) :
		"""SNVsFilter is a fct tha takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		
		if symbol not in self.genes.keys() :
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
		
	#def loadGene_bck(self, symbolOrId) :
	#	"""Loads a gene and returns it"""
	#	if symbolOrId not in self.genes.keys() :
	#
	#		try :
	#			f = open('ncbi/%s/sequences/genes/chr%s_%s.gtf'%(self.genome.name, self.number, symbolOrId))
	#			geneData = f.read()
	#			f.close()
	#		except :
	#			geneData = "".join(re.findall('.+"%s".+\n'%(symbolOrId), self.gtf))
	#			f = open('ncbi/%s/sequences/genes/chr%s_%s.gtf'%(self.genome.name, self.number, symbolOrId), 'w')
	#			f.write(geneData)
	#			f.close()
	#		#print geneData
	#		self.genes[symbolOrId] = Gene(symbolOrId, self, geneData)
	#	
	#	return self.genes[symbolOrId]
	
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
				#print '----', snps
				data = None
				if snps != None :
					for snp in snps:
						#print '-----'
						#print start, end, snp
						if fct(snp) :
							if data == None :
								data = list(self.data[start:end])
							pos = snp['pos'] - start#-1
							#print start, end, snp['pos'], start, pos, len(data)
							data[pos] = snp['max_gt']
							#print snp['max_gt']
				if data != None :
					return ''.join(data)
				else :
					return self.data[start:end]
	
	def getPolymorphismsInRange(self, x1, x2) :
		return self.casavaSNPs.findSnpsInRange(x1, x2)
		
	def getSequence_bck(self, x1, x2, SNVsFilter = None) :
		"""SNVsFilter is a fct tha takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		
		if x1 > x2 :
			start, end = x2, x1 
		else :
			start, end = x1, x2
			
		#print "str:", self.data[end:start], len(self.data), start, end, self.data[end:start]
		if not self.isLight :
			if start <= end :
				return self.data[start:end]
			return self.data[end:start]
		else :
			if (SNVsFilter != None) :
				fct = SNVsFilter
			else :
				fct = defaultSNVsFilter
			
			snps = self.casavaSNPs.findSnpsInRange(start, end)
			#print snps
			data = list(self.data[start:end])
			data2 = list(self.data[start:end])
			if snps != None :
				for snp in snps:					
					if fct(snp) :
						pos = snp['pos'] - start-1
						data[pos] = snp['max_gt']
					if defaultSNVsFilter(snp):
						pos = snp['pos'] - start-1
						data2[pos] = snp['max_gt']
						
			if data != data2 :
				for i in range(len(data)) :
					if data[i] != data2[i] :
						print i, data[i], data2[i]
					
			if data != None :
				return ''.join(data)
			else :
				return self.data[start:end]
				
	def getRandomCDS(self) :
		cds = CDS()
		found = False
		while not found :
			#print l, len(self.gtfLines)-1
			#print self.gtfLines[l]
			found = cds.setFromGTFLine(random.sample(self.gtfLines, 1)[0])
		
		return cds

	def getRegionIndex(self, name = '') :
		return None
		#return RegionIndex(self, name)
	
	def getSNPsInRange(self, x1, x2) :
		"""DB snps"""
		snps =[]
		line = self.SNPIndexTree.getStartLine(x1+1)#je compte a 0 ncbi a 1
		snp = self.materialiseSNP(line)
		while snp.x1 < x2 : #<=
			if x1 <= snp.x1 :
				snps.append(snp)
			line += 1
			snp = self.materialiseSNP(line)
		
		return snps
	
	def getSNVsInRange(self, x1, x2, SNVsFilter = None) :
		"""SNVsFilter is a fct that takes a CasavaSnp as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		res = []
		if x1 > x2 :
			start, end = x2, x1 
		else :
			start, end = x1, x2
			
		if not self.isLight :
			return res
		else :
			if (SNVsFilter != None) :
				fct = SNVsFilter
			else :
				fct = defaultSNVsFilter
			
			snps = self.casavaSNPs.findSnpsInRange(start, end)
			
			if snps != None :
				for snp in snps:
					if fct(snp) :
						res.append(snp)
		return res
		
	#def translateSequence(self, start, end) :
	#	return tools.translateDNA(self.getSequence(start, end))
	
	#def translateSequence_6Read(self, start, end) :
	#	return tools.translateDNA_6Read(self.getSequence(start, end))
		
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
		#"""Returns a gene"""
		return self.genes[i]
		#return self.data[i]
		
	def __len__(self):
		return len(self.data)
