import configuration as conf

from pyGenoObjectBases import *
from SNP import SNP_INDEL

import rabaDB.fields as rf

class Gene_Raba(pyGenoRabaObject) :
	"""The wrapped Raba object that really holds the data"""
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	id = rf.Primitive()
	name = rf.Primitive()
	strand = rf.Primitive()
	biotype = rf.Primitive()
	
	start = rf.Primitive()
	end = rf.Primitive()
	
	genome = rf.RabaObject('Genome_Raba')
	chromosome = rf.RabaObject('Chromosome_Raba')

	def _curate(self) :
		self.name = self.name.upper()

class Gene(pyGenoRabaObjectWrapper) :
	"""The wrapper for playing with genes"""
	
	_wrapped_class = Gene_Raba

	def __init__(self, *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)

	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		if issubclass(objectType, SNP_INDEL) :
			f = RabaQuery(objectType, namespace = self._wrapped_class._raba_namespace)
			coolArgs['species'] = self.genome.species
			coolArgs['chromosomeNumber'] = self.chromosome.number
			coolArgs['start >='] = self.start
			coolArgs['start <'] = self.end
		
			if len(args) > 0 and type(args[0]) is types.ListType :
				for a in args[0] :
					if type(a) is types.DictType :
						f.addFilter(**a)
			else :
				f.addFilter(*args, **coolArgs)

			return f
		
		return pyGenoRabaObjectWrapper._makeLoadQuery(self, objectType, *args, **coolArgs)
	
	def __str__(self) :
		return "Gene, name: %s, id: %s, strand: '%s' > %s" %(self.name, self.id, self.strand, str(self.chromosome))
