import types
import configuration as conf
import pyGeno.tools.UsefulFunctions as uf
from pyGenoObjectBases import *

from pyGeno.Chromosome import Chromosome
from pyGeno.Gene import Gene
from pyGeno.Transcript import Transcript
from pyGeno.Protein import Protein
from pyGeno.Exon import Exon
from SNP import *

import rabaDB.fields as rf

class Genome_Raba(pyGenoRabaObject) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	#_raba_not_a_singleton = True #you can have several instances of the same genome but they all share the same location in the database

	name = rf.Primitive()
	specie = rf.Primitive()

	source = rf.Primitive()
	packageInfos = rf.Primitive()

	def _curate(self) :
		self.specie = self.specie.lower()

	def getSequencePath(self) :
		return conf.getGenomeSequencePath(self.specie, self.name)

	def getReferenceSequencePath(self) :
		return conf.getReferenceGenomeSequencePath(self.specie)

	def __len__(self) :
		"""Size of the genome in pb"""
		l = 0
		for c in self.chromosomes :
			l +=  len(c)

		return l

class Genome(pyGenoRabaObjectWrapper) :
	_wrapped_class = Genome_Raba

	def __init__(self, SNPs = None, SNPFilter = defaultSNPFilter,  *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)

		if type(SNPs) is types.StringType :
			self.SNPsSets = [SNPs]
		else :
			self.SNPsSets = SNPs
		
		self.SNPFilter = SNPFilter
		self.SNPTypes = {}
		
		if SNPs is not None :
			f = RabaQuery(SNPMaster, namespace = self._raba_namespace)
			for se in self.SNPsSets :
				f.addFilter(setName = se, specie = self.specie)

			res = f.run()
			if res is None or len(res) < 1 :
				raise ValueError("There's no set of SNPs that goes by the name of %s for specie %s" % (SNPs, self.specie))

			for s in res :
				self.SNPTypes[s.setName] = s.SNPType

	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		if issubclass(objectType, SNP_INDEL) :
			f = RabaQuery(objectType, namespace = self._wrapped_class._raba_namespace)
			coolArgs['specie'] = self.specie

			if len(args) > 0 and type(args[0]) is types.ListType :
				for a in args[0] :
					if type(a) is types.DictType :
						f.addFilter(**a)
			else :
				f.addFilter(*args, **coolArgs)

			return f
		
		return pyGenoRabaObjectWrapper._makeLoadQuery(self, objectType, *args, **coolArgs)
	
	def __str__(self) :
		return "Genome: %s/%s SNPs: %s" %(self.specie, self.name, self.SNPTypes)
