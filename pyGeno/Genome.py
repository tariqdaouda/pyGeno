import types
import configuration as conf
import pyGeno.tools.UsefulFunctions as uf
from pyGenoObjectBases import *

from Chromosome import Chromosome
from Gene import Gene
from Transcript import Transcript
from Protein import Protein
from Exon import Exon
import SNPFiltering as SF
from SNP import *

import rabaDB.fields as rf

def getGenomeList() :
	"""Return the names of all imported genomes"""
	import rabaDB.filters as rfilt
	f = rfilt.RabaQuery(Genome_Raba)
	names = []
	for g in f.iterRun() :
		names.append(g.name)
	return names

class Genome_Raba(pyGenoRabaObject) :
	"""The wrapped Raba object that really holds the data"""
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	#_raba_not_a_singleton = True #you can have several instances of the same genome but they all share the same location in the database

	name = rf.Primitive()
	species = rf.Primitive()

	source = rf.Primitive()
	packageInfos = rf.Primitive()

	def _curate(self) :
		self.species = self.species.lower()

	def getSequencePath(self) :
		return conf.getGenomeSequencePath(self.species, self.name)

	def getReferenceSequencePath(self) :
		return conf.getReferenceGenomeSequencePath(self.species)

	def __len__(self) :
		"""Size of the genome in pb"""
		l = 0
		for c in self.chromosomes :
			l +=  len(c)

		return l

class Genome(pyGenoRabaObjectWrapper) :	
	"""
	This is the entry point to pyGeno::
		
		myGeno = Genome(name = 'GRCh37.75', SNPs = ['RNA_S1', 'DNA_S1'], SNPFilter = MyFilter)
		for prot in myGeno.get(Protein) :
			print prot.sequence
	
	"""
	_wrapped_class = Genome_Raba

	def __init__(self, SNPs = None, SNPFilter = None,  *args, **kwargs) :
		
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)

		if type(SNPs) is types.StringType :
			self.SNPsSets = [SNPs]
		else :
			self.SNPsSets = SNPs
		
		if SNPFilter is None :
			self.SNPFilter = SF.DefaultSNPFilter()
		else :
			if issubclass(SNPFilter.__class__, SF.SNPFilter) :
				self.SNPFilter = SNPFilter
			else :
				raise ValueError("The value of 'SNPFilter' is not an object deriving from a subclass of SNPFiltering.SNPFilter. Got: '%s'" % SNPFilter)

		self.SNPTypes = {}
		
		if SNPs is not None :
			f = RabaQuery(SNPMaster, namespace = self._raba_namespace)
			for se in self.SNPsSets :
				f.addFilter(setName = se, species = self.species)

			res = f.run()
			if res is None or len(res) < 1 :
				raise ValueError("There's no set of SNPs that goes by the name of %s for species %s" % (SNPs, self.species))

			for s in res :
				self.SNPTypes[s.setName] = s.SNPType

	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		if issubclass(objectType, SNP_INDEL) :
			# conf.db.enableDebug(True)
			f = RabaQuery(objectType, namespace = self._wrapped_class._raba_namespace)
			coolArgs['species'] = self.species

			if len(args) > 0 and type(args[0]) is types.ListType :
				for a in args[0] :
					if type(a) is types.DictType :
						f.addFilter(**a)
			else :
				f.addFilter(*args, **coolArgs)

			return f
		
		return pyGenoRabaObjectWrapper._makeLoadQuery(self, objectType, *args, **coolArgs)
	
	def __str__(self) :
		return "Genome: %s/%s SNPs: %s" %(self.species, self.name, self.SNPTypes)
