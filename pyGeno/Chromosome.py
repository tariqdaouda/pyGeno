#import copy
#import types
#from tools import UsefulFunctions as uf

from types import *
import configuration as conf
from pyGenoObjectBases import *

from SNP import *
import SNPFiltering as SF

from rabaDB.filters import RabaQuery
import rabaDB.fields as rf

from tools.SecureMmap import SecureMmap as SecureMmap
from tools import SingletonManager

import pyGeno.configuration as conf

class ChrosomeSequence(object) :
	"""Represents a chromosome sequence. If 'refOnly' no ploymorphisms are applied and the ref sequence is always returned"""

	def __init__(self, data, chromosome, refOnly = False) :
		
		self.data = data
		self.refOnly = refOnly
		self.chromosome = chromosome
		self.setSNPFilter(self.chromosome.genome.SNPFilter)
	
	def setSNPFilter(self, SNPFilter) :
		self.SNPFilter = SNPFilter
	
	def getSequenceData(self, slic) :
		data = self.data[slic]
		SNPTypes = self.chromosome.genome.SNPTypes
		if SNPTypes is None or self.refOnly :
			return data
		
		iterators = []
		for setName, SNPType in SNPTypes.iteritems() :
			f = RabaQuery(str(SNPType), namespace = self.chromosome._raba_namespace)
			
			chromosomeNumber = self.chromosome.number

			if chromosomeNumber == 'MT':
				chromosomeNumber = 'M'
			
			f.addFilter({'start >=' : slic.start, 'start <' : slic.stop, 'setName' : str(setName), 'chromosomeNumber' : chromosomeNumber})
			# conf.db.enableDebug(True)
			iterators.append(f.iterRun(sqlTail = 'ORDER BY start'))
		
		if len(iterators) < 1 :
			return data
		
		polys = {}
		for iterator in iterators :
			for poly in iterator :
				if poly.start not in polys :
					polys[poly.start] = {poly.setName : poly}
				else :
					try :
						polys[poly.start][poly.setName].append(poly)
					except :
						polys[poly.start][poly.setName] = [polys[poly.start][poly.setName]]
						polys[poly.start][poly.setName].append(poly)
						
		data = list(data)
		for start, setPolys in polys.iteritems() :
			
			seqPos = start - slic.start
			sequenceModifier = self.SNPFilter.filter(self.chromosome, **setPolys)
			# print sequenceModifier.alleles
			if sequenceModifier is not None :
				if sequenceModifier.__class__ is SF.SequenceDel :
					seqPos = seqPos + sequenceModifier.offset
					#To avoid to change the length of the sequence who can create some bug or side effect
					data[seqPos:(seqPos + sequenceModifier.length)] = [''] * sequenceModifier.length
				elif sequenceModifier.__class__ is SF.SequenceSNP :
					data[seqPos] = sequenceModifier.alleles
				elif sequenceModifier.__class__ is SF.SequenceInsert :
					seqPos = seqPos + sequenceModifier.offset
					data[seqPos] = "%s%s" % (data[seqPos], sequenceModifier.bases)
				else :
					raise TypeError("sequenceModifier on chromosome: %s starting at: %s is of unknown type: %s" % (self.chromosome.number, snp.start, sequenceModifier.__class__))

		return data
	
	def _getSequence(self, slic) :
		return ''.join(self.getSequenceData(slic))

	def __getitem__(self, i) :
		return self._getSequence(i)

	def __len__(self) :
		return self.chromosome.length

class Chromosome_Raba(pyGenoRabaObject) :
	"""The wrapped Raba object that really holds the data"""
	
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
	"""The wrapper for playing with Chromosomes"""
	
	_wrapped_class = Chromosome_Raba

	def __init__(self, *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)

		path = '%s/chromosome%s.dat'%(self.genome.getSequencePath(), self.number)
		if not SingletonManager.contains(path) :
			datMap = SingletonManager.add(SecureMmap(path), path)
		else :
			datMap = SingletonManager.get(path)
			
		self.sequence = ChrosomeSequence(datMap, self)
		self.refSequence = ChrosomeSequence(datMap, self, refOnly = True)
		self.loadSequences = False

	def getSequenceData(self, slic) :
		return self.sequence.getSequenceData(slic)

	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		if issubclass(objectType, SNP_INDEL) :
			f = RabaQuery(objectType, namespace = self._wrapped_class._raba_namespace)
			coolArgs['species'] = self.genome.species
			coolArgs['chromosomeNumber'] = self.number

			if len(args) > 0 and type(args[0]) is types.ListType :
				for a in args[0] :
					if type(a) is types.DictType :
						f.addFilter(**a)
			else :
				f.addFilter(*args, **coolArgs)

			return f
		
		return pyGenoRabaObjectWrapper._makeLoadQuery(self, objectType, *args, **coolArgs)

	def __getitem__(self, i) :
		return self.sequence[i]

	def __str__(self) :
		return "Chromosome: number %s > %s" %(self.wrapped_object.number, str(self.wrapped_object.genome))
