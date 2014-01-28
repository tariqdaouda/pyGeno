from types import *
import configuration as conf
from pyGenoObject import *

#from rabaDB.setup import *
#RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
#from rabaDB.Raba import *
#from rabaDB.filters import RabaQuery
import rabaDB.fields as rf

#from Gene import Gene
#from SNP import *

#from tools.SegmentTree import SegmentTree as SegmentTree
from tools.SecureMmap import SecureMmap as SecureMmap
from tools import UsefulFunctions as uf
from tools import SingletonManager
import types

def defaultSNPsFilter(casavaSnp) :
	"""The default rule to decide wether to take the most probable genotype or the
	reference, always returns true"""
	return True

class ChrosomeSequence(object) :

	def __init__(self, data, chromosome) :
		self.data = data
		self.chromosome = chromosome
		self.SNPsFilter = defaultSNPsFilter

	def setSNPFilter(self, SNPsFilter) :
		self.SNPsFilter = SNPsFilter

	def _getSequence(self, slic) :
		assert type(slic) is SliceType

		if self.chromosome.dataType == 'heavy' :
			return self.data[slic]
		else :
			snps = self.chromosome.get(self.chromosome.dataType, {'pos >=' : slic[0], 'pos <' : slic[1]})
			data = self.data[slic]
			if len(snps) < 1 :
				return data
			for snp in snps:
				if self.SNPsFilter(snp) :
					if type(data) is not ListType :
						data = list(data)
					posSeq = snp['pos'] - start#-1
					snp['max_gt'] = uf.getPolymorphicNucleotide(snp['max_gt'])
					data[posSeq] = snp['max_gt']

			if type(data) is ListType :
				return ''.join(data)
			else :
				return data

	def __getitem__(self, tralala) :
		return self._getSequence(tralala)

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
	dataType = rf.Primitive() #'flat' => for dat files on drive, or name of the polymoprhism's rabaclass : ex 'CasavaSNP'

	genome = rf.RabaObject('Genome')
	#genes = rf.Relation('Gene')

	#_raba_uniques = [('genome', 'number')]

	def __init__(self, importing = False, SNPsFilter = defaultSNPsFilter) :
		"""SNVsFilter is a fct that takes a SNP as input a returns true if it correpsond to the rule.
		If left to none Chromosome.defaulSNVsFilter is used. This parameter has no effect if the genome is not light
		(contains the sequences for all chros)"""
		if self.number != None :
			self.number = str(self.number)

		if not importing :
			if self.dataType == 'heavy' :
				path = '%s/chromosome%s.dat'%(self.genome.getSequencePath(), self.number)
				self.sequence = ChrosomeSequence(SingletonManager.add(SecureMmap(path), path), self)

			self.SNPsFilter = SNPsFilter

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

	def __getitem__(self, i) :
		return self.sequence[i]

	def __len__(self) :
		return self.length

	def __str__(self) :
		return "Chromosome: number %s > %s" %(self.number, str(self.genome))
