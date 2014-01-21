import os
import configuration as conf
from pyGenoObject import *

#from rabaDB.setup import *
#RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
#from rabaDB.Raba import *
#from rabaDB.filters import RabaQuery
import rabaDB.fields as rf

class Genome(pyGenoObject) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	name = rf.Primitive()
	specie = rf.Primitive()
	#chromosomes = rf.Relation('Chromosome')
	snpsType = rf.Primitive()

	source = rf.Primitive()
	packageInfos = rf.Primitive()

	_raba_uniques = [('name', 'specie')]

	def __init__(self) :
		pass

	def _curate(self) :
		pass

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

	def __str__(self) :
		return "Genome: %s/%s" %(self.specie, self.name)

