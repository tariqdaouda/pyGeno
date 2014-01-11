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

def defaultSNVsFilter(casavaSnp) :
	"""The default rule to decide wether to take the most probable genotype or the
	reference, always returns true"""
	return True

def defaultDbSNPsFilter(dbSnp) :
	"""The default rule to decide wether to take the most probable genotype or the
	snp, always returns true"""
	return True
	
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
	genes = rf.Relation('Gene')
	
	_raba_uniques = [('genome', 'number')]
	
	def __init__(self) :
		if self.number != None :
			self.number = str(self.number)
			self._loadSequence()
			if self.dataType == 'flat' :
				path = '%s/chr%s.dat'%(self.genome.getSequencePath(), self.number)
				self.data = SingletonManager.add(SecureMmap(path), path)

	def save(self) :
		if  self.x2 != None and self.x1 != None :
			self.length = self.x2-self.x1
		if self.number != None :
			self.number = self.number.upper()
		
		pyGenoObject.save(self)
	
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
		
		assert type(x1) is IntType
		assert type(x2) is IntType
		
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
				data = None
				if snps != None :
					for snp in snps:
						if fct(snp) :
							if data == None :
								data = list(self.data[start:end])
							pos = snp['pos'] - start#-1
							snp['max_gt'] = uf.getPolymorphicNucleotide(snp['max_gt'])
							data[pos] = snp['max_gt']

				if data != None :
					return ''.join(data)
				else :
					return self.data[start:end]
	
	def getPolymorphismsInRange(self, x1, x2) :
		return self.casavaSNPs.findSnpsInRange(x1, x2)
	
	def stringFind(self, sequence) :
		return self.data.find(sequence)

	def pluck(self) :
		"""Returns a plucked object. Plucks the chromosome off the tree, set the value of self.genome into str(self.genome). This effectively disconnects the object and
		makes it much more lighter in case you'd like to pickle it"""
		e = copy.copy(self)
		e.genome = str(self.genome)
		return e
		
	def __getitem__(self, i) :
		return self.genes[i]
		
	def __len__(self) :
		return self.length

	def __str__(self) :
		return "Chromosome: number %s > %s" %(self.number, str(self.genome))
