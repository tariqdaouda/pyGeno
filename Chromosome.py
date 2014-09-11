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

	def __init__(self, data, chromosome, refOnly = False) :
		"""represents a chromosome sequence. if refOnly no ploymorphisms are applied and the ref sequence is always returned"""

		self.data = data
		self.refOnly = refOnly
		self.chromosome = chromosome
		self.SNPFilter = self.chromosome.genome.SNPFilter
	
	def setSNPFilter(self, SNPFilter) :
		self.SNPFilter = SNPFilter
	
	def _getSequence(self, slic) :
		#~ print slic
		#~ assert type(slic) is SliceType
		#~ 
		data = self.data[slic]
		SNPTypes = self.chromosome.genome.SNPTypes
		
		if SNPTypes is None or self.refOnly :
			return data
		
		iterators = []
		for setName, SNPType in SNPTypes.iteritems() :
			f = RabaQuery(str(SNPType), namespace = self.chromosome._raba_namespace)
			f.addFilter({'start >=' : slic.start, 'start <' : slic.stop, 'setName' : str(setName), 'chromosomeNumber' : self.chromosome.number})
			#conf.db.enableDebug(True)
			iterators.append(f.iterRun(sqlTail = 'ORDER BY start'))
		
		if len(iterators) < 1 :
			return data
		
		polys = {}
		for iterator in iterators :
			for poly in iterator : 
				if poly.start not in polys :
					polys[poly.start] = {poly.setName : poly}
				else :
					polys[poly.start][poly.setName] = poly
		
		data = list(data)
		for start, setPolys in polys.iteritems() :
			
			seqPos = start - slic.start
			sequenceSNP = self.SNPFilter(chromosome = self.chromosome, **setPolys)
			
			if sequenceSNP.length < 1 :
				raise TypeError("SequenceSNP of chromosome: %s starting at: %s has a .length < 1 (%s)" % (self.chromosome.number, start, sequenceSNP.length))
			
			if sequenceSNP.alleles not in uf.polymorphicNucleotides and sequenceSNP.alleles not in uf.nucleotides :
				raise TypeError("SequenceSNP of chromosome: %s starting at: %s has invalid alleles" % (self.chromosome.number, start, sequenceSNP.alleles))
			
			if sequenceSNP.type is SequenceSNP_INDEL.DeletionType :
				data = data[:seqPos-1] + data[seqPos-1 + sequenceSNP.length] 
			elif sequenceSNP.type is SequenceSNP_INDEL.InsertionType or sequenceSNP.type is SequenceSNP_INDEL.SNPType :
				data[seqPos] = sequenceSNP.alleles
			else :
				raise TypeError("SequenceSNP of chromosome: %s starting at: %s is of unknown type: %s" % (self.chromosome.number, snp.start, sequenceSNP.type))
		
		return ''.join(data)

	def __getitem__(self, i) :
		return self._getSequence(i)
		#return self.data[i]

	def __len__(self) :
		return self.chromosome.length

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
		self.refSequence = ChrosomeSequence(SingletonManager.add(SecureMmap(path), path), self, refOnly = True)
		self.loadSequences = False
	
	def stringFind(self, sequence) :
		return self.sequence.find(sequence)

	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		if issubclass(objectType, SNP_INDEL) :
			f = RabaQuery(objectType, namespace = self._wrapped_class._raba_namespace)
			coolArgs['chromosomeNumber'] = self.number

			if len(args) > 0 and type(args[0]) is types.ListType :
				for a in args[0] :
					if type(a) is types.DictType :
						f.addFilter(**a)
			else :
				f.addFilter(*args, **coolArgs)

			return f
		
		return pyGenoRabaObjectWrapper._makeLoadQuery(self, objectType, *args, **coolArgs)

	def __str__(self) :
		return "Chromosome: number %s > %s" %(self.wrapped_object.number, str(self.wrapped_object.genome))
