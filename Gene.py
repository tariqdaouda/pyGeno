import configuration as conf
from pyGenoObject import *

#from rabaDB.setup import *
#RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
#from rabaDB.Raba import *
#from rabaDB.filters import RabaQuery
import rabaDB.fields as rf

#from tools import UsefulFunctions as uf
#from Protein import Protein
#from Exon import Exon
#from Transcript import Transcript

#from tools.BinarySequence import NucBinarySequence

class Gene(pyGenoObject) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	
	id = rf.Primitive()
	name = rf.Primitive()
	strand = rf.Primitive()
	biotype = rf.Primitive()
	
	genome = rf.RabaObject('Genome')
	chromosome = rf.RabaObject('Chromosome')
	transcripts = rf.Relation('Transcript')
	exons = rf.Relation('Exon')
	
	_raba_uniques = [('genome', 'id')]
	
	def _curate(self) :
		self.name = self.name.upper()
	
	def __init__(self, *args, **fieldsSet) :
		pass
	
	def pluck(self) :
		"""Returns a plucked object. Plucks the gene off the tree, set the value of self.chromosome into str(self.chromosome). This effectively disconnects the object and
		makes it much more lighter in case you'd like to pickle it"""
		e = copy.copy(self)
		e.chromosome = str(self.chromosome)
		return e
	
	def __str__(self) :
		return "Gene, name: %s, id: %s, strand: '%s' > %s" %(self.name, self.id, self.strand, str(self.chromosome))
