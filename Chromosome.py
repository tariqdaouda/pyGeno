import copy
from types import *
import configuration as conf
from pyGenoObject import *

#from rabaDB.setup import *
#RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
#from rabaDB.Raba import *
from rabaDB.filters import RabaQuery
import rabaDB.fields as rf

#from Gene import Gene
#from SNP import *

#from tools.SegmentTree import SegmentTree as SegmentTree
from tools.SecureMmap import SecureMmap as SecureMmap
from tools import UsefulFunctions as uf
from tools import SingletonManager
import types


class ChrosomeSequence(object) :

	def __init__(self, data, chromosome) :
		self.data = data
		self.chromosome = chromosome
		self.SNPsFilter = self.chromosome.genome.SNPsFilter

	def setSNPFilter(self, SNPsFilter) :
		self.SNPsFilter = SNPsFilter

	def _getSequence(self, slic) :
		"returns a sequence including SNPs as filtered by the genome SNPFilter function"
		assert type(slic) is SliceType

		data = self.data[slic]
		SNPTypes = self.chromosome.genome.SNPTypes
		if SNPTypes != None :
			resSNPs = []
			for SNPType in SNPTypes.itervalues() :
				f = RabaQuery(snpType, namespace = self.chromosome._raba_namespace)
				f.addFilter({'start >=' : slic[0], 'start < ' : slic[1], 'setName' : se, 'chromosomeNumber' : self.chromosome.number})
				resSNPs.append(f.run(sqlTail = 'ORDER BY position'))

			if len(resSNPs) == 1 :
				for SNP in resSNPs[0] :
					filtSNP = self.SNPsFilter(SNP)
					if filtSNP != None :
						if type(data) is not ListType :
							data = list(data)
						posSeq = filtSNP.start - slic[0]#-1
						filtSNP.alleles = uf.getPolymorphicNucleotide(filtSNP.alleles)
						data[posSeq] = filtSNP.alleles
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

	def _mixSNPs(*snpsSets) :
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

class Chromosome(pyGenoObject) :
	"""A class that represents a Chromosome
	Attention: private region support en retard par rapport au public"""
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	header = rf.Primitive()
	number = rf.Primitive()
	#x1, x2 are the prosition of the chromosome in the genome
	x1 = rf.Primitive()
	x2 = rf.Primitive()
	length = rf.Primitive()
	#dataType = rf.Primitive() #'flat' => for dat files on drive, or name of the polymoprhism's rabaclass : ex 'CasavaSNP'

	genome = rf.RabaObject('Genome')
	#genes = rf.Relation('Gene')

	#_raba_uniques = [('genome', 'number')]

	def __init__(self, importing = False) :
		"""SNVsFilter is a fct that takes a SNP as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		if self.number != None :
			self.number = str(self.number)

		if not importing :
			path = '%s/chromosome%s.dat'%(self.genome.getSequencePath(), self.number)
			self.sequence = ChrosomeSequence(SingletonManager.add(SecureMmap(path), path), self)

		#	self.SNPsFilter = SNPsFilter

	def _curate(self) :
		if  self.x2 != None and self.x1 != None :
			self.length = self.x2-self.x1
		if self.number != None :
			self.number =  str(self.number).upper()

	def setSNPFilter(self, SNPsFilter) :
		self.SNPsFilter = SNPsFilter

	def getPolymorphismsInRange(self, x1, x2) :
		return self.casavaSNPs.findSnpsInRange(x1, x2)

	def stringFind(self, sequence) :
		return self.sequence.find(sequence)

	def pluck(self) :
		"""Returns a plucked object. Plucks the chromosome off the tree, set the value of self.genome into str(self.genome). This effectively disconnects the object and
		makes it much more lighter in case you'd like to pickle it"""
		e = copy.copy(self)
		e.genome = str(self.genome)
		return e

	def __len__(self) :
		return self.length

	def __str__(self) :
		return "Chromosome: number %s > %s" %(self.number, str(self.genome))
