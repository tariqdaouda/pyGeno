import configuration as conf
from pyGenoObjectBases import *

import rabaDB.fields as rf

class Gene_Raba(pyGenoRabaObject) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	id = rf.Primitive()
	name = rf.Primitive()
	strand = rf.Primitive()
	biotype = rf.Primitive()

	genome = rf.RabaObject('Genome_Raba')
	chromosome = rf.RabaObject('Chromosome_Raba')

	def _curate(self) :
		self.name = self.name.upper()

class Gene(pyGenoRabaObjectWrapper) :
	_wrapped_class = Gene_Raba

	def __init__(self, *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)

	def __str__(self) :
		return "Gene, name: %s, id: %s, strand: '%s' > %s" %(self.name, self.id, self.strand, str(self.chromosome))
