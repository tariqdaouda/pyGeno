import copy
from types import *
import configuration as conf
from pyGenoObjectBases import *

from SNP import *
from rabaDB.filters import RabaQuery
import rabaDB.fields as rf

from tools.SecureMmap import SecureMmap as SecureMmap
from tools import UsefulFunctions as uf
from tools import SingletonManager
import types

import pyGeno.configuration as conf

class ChrosomeSequence(object) :

	def __init__(self, data, chromosome) :
		self.data = data
		self.chromosome = chromosome
		self.SNPsFilter = self.chromosome.genome.SNPsFilter
	
	def setSNPFilter(self, SNPsFilter) :
		self.SNPsFilter = SNPsFilter

	def _getSequence(self, slic) :
		assert type(slic) is SliceType
		
		data = self.data[slic]
		SNPTypes = self.chromosome.genome.SNPTypes
		
		if SNPTypes is not None :
			iterators = {}
			
			for setName, SNPType in SNPTypes.iteritems() :
				f = RabaQuery(str(SNPType), namespace = self.chromosome._raba_namespace)
				f.addFilter({'start >=' : slic.start, 'start <' : slic.stop, 'setName' : str(setName), 'chromosomeNumber' : self.chromosome.number})
				#conf.db.enableDebug(True)
				iterators[setName] = f.iterRun(sqlTail = 'ORDER BY start')
			
			if len(iterators) < 1 :
				return data
			
			data = list(data)
			
			lowStack = {}
			lowStackSortedStarts = []
			setsInLowStack = set()
			deletions = [] #tuples (seqPos, len)
			
			while len(iterators) > 0 :
				
				itsWithoutSNPs = []
				for setName, iterator in iterators.iteritems() :
					try :
						if setName not in setsInLowStack :
							snp = iterator.next()
							if snp.start not in lowStack :
								lowStack[snp.start] = {}
								lowStackSortedStarts.append(snp.start)
							lowStack[snp.start][setName] = snp
							setsInLowStack.add(setName)
					except StopIteration :
						itsWithoutSNPs.append(setName)
				
				if len(lowStackSortedStarts) < 1 :
					break # no snps in that range
				
				lowStackSortedStarts.sort()
				pos = lowStackSortedStarts.pop()
				posSeq = snp.start - slic.start#-1
				
				#print snp.alt
				refAllele = data[posSeq:posSeq +len(snp.alt)]
				alleles = self.SNPsFilter(refAllele = refAllele, **lowStack[pos])
				for setName in lowStack[pos] :
					setsInLowStack.remove(setName)
				del(lowStack[pos])
				if len(refAllele) > len(alleles) :
					deletions.append(posSeq, len(refAllele) - 1) #-1 to keep the inserted allele
				#print posSeq, posSeq, data[posSeq], alleles
				data[posSeq] = str(alleles)
				
				for setName in itsWithoutSNPs :
					del(iterators[setName])
				
		for posSeq, length in deletions :
			for i in range(posSeq + 1, length) :
				del(data[i])
		
		if type(data) is ListType :
			data = ''.join(data)
			#print data
		
		return data

	def _getSequence_bck(self, slic) :
		"returns a sequence including SNPs as filtered by the genome SNPFilter function"
		assert type(slic) is SliceType

		data = self.data[slic]
		SNPTypes = self.chromosome.genome.SNPTypes

		if SNPTypes != None :
			resSNPs = []
			for setName, SNPType in SNPTypes.iteritems() :
				f = RabaQuery(str(SNPType), namespace = self.chromosome._raba_namespace)
				f.addFilter({'start >=' : slic.start, 'start <' : slic.stop, 'setName' : str(setName), 'chromosomeNumber' : self.chromosome.number})
				resSNPs.append(f.run(sqlTail = 'ORDER BY start'))

			if len(resSNPs) == 1 :
				for SNP in resSNPs[0] :
					filtSNP = self.SNPsFilter(SNP)
					if filtSNP != None :
						if type(data) is not ListType :
							data = list(data)
						posSeq = filtSNP.start - slic.start#-1
						filtSNP.alleles = uf.getPolymorphicNucleotide(filtSNP.alt)
						data[posSeq] = str(filtSNP.alt)
						
			elif len(resSNPs) > 1 :
				for SNP in self._mixSNPs(*resSNPs) :
					filtSNP = self.SNPsFilter(**SNP)
					if filtSNP != None :
						if type(data) is not ListType :
							data = list(data)
						posSeq = filtSNP.start - slic[0]#-1
						filtSNP.alleles = uf.getPolymorphicNucleotide(filtSNP.alleles)
						data[posSeq] = filtSNP.alleles

			if type(data) is ListType :
				return ''.join(data)

		return data

	def _mixSNPs_bck(*snpsSets) :
		"""takes several snp sets and return an iterator of values {setName1 : snp, setName2 : snp, setName3 : None, ...}.
		the dict a intended for the SNPFilter function, each one correponds to a single position in the chromosome"""

		positions = {}
		empty = {}
		for snpsSet in snpsSets :
			empty[snpsSet[0].setName] = None

		for snpsSet in snpsSets :
			for snp in snpsSet :
				if snp.start not in position :
					positions[snp.start] = copy.copy(empty)
				positions[snp.start][empty[snp.setName]] = snp

		for v in positions.itervalues() :
			yield v

	def __getitem__(self, i) :
		return self._getSequence(i)

class Chromosome_Raba(pyGenoRabaObject) :
	"""A class that represents a persistent Chromosome"""
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	header = rf.Primitive()
	number = rf.Primitive()
	start = rf.Primitive()
	end = rf.Primitive()
	length = rf.Primitive()

	genome = rf.RabaObject('Genome_Raba')

	def _curate(self) :
		if  self.end != None and self.start != None :
			self.length = self.end-self.start
		if self.number != None :
			self.number =  str(self.number).upper()

class Chromosome(pyGenoRabaObjectWrapper) :

	_wrapped_class = Chromosome_Raba

	def __init__(self, *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)

		path = '%s/chromosome%s.dat'%(self.genome.getSequencePath(), self.number)
		self.sequence = ChrosomeSequence(SingletonManager.add(SecureMmap(path), path), self)
		self.loadSequences = False
	
	def stringFind(self, sequence) :
		return self.sequence.find(sequence)

	def __str__(self) :
		return "Chromosome: number %s > %s" %(self.wrapped_object.number, str(self.wrapped_object.genome))
