import os
import configuration as conf
import pyGeno.tools.UsefulFunctions as uf
from pyGenoObjectBases import *

from pyGeno.Chromosome import Chromosome
from pyGeno.Gene import Gene
from pyGeno.Transcript import Transcript
from pyGeno.Protein import Protein
from pyGeno.Exon import Exon
from SNP import *

#from rabaDB.setup import *
#RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
#from rabaDB.Raba import *
#from rabaDB.filters import RabaQuery
import rabaDB.fields as rf

"""
class SequenceSNP_INDEL(object) :
	"This SNP is intended to be embeded in a sequence and has no persistency. It's alleles can be a mix of alleles from SNPs of different sets"
	def __init__(self, alleles, start, end) :
		self.alleles = alleles
		self.start = start
		self.end = end
		self.sourceSNPs = []

	def addSourceSNP(self, sourceSNP) :
		"persistent SNPs whose alleles have been mixed to make self.alleles"
		self.sourceSNPs.append(sourceSNP)
"""

def defaultSNPsFilter(refAllele, *args, **kwargs) :
	"Default function for filtering snp, simply does not filter anything. Doesn't work with indels"
	alleles = []
	for snp in args :
		alleles.append(snp.alt)
	
	for k, snp in kwargs.iteritems() :
		alleles.append(snp.alt)
	
	return uf.encodePolymorphicNucleotide(''.join(alleles))

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

	def __init__(self, SNPs = None, SNPsFilter = defaultSNPsFilter,  *args, **kwargs) :
		pyGenoRabaObjectWrapper.__init__(self, *args, **kwargs)

		self.SNPsSets = SNPs
		self.SNPsFilter = SNPsFilter
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

	def __str__(self) :
		return "Genome: %s/%s SNPs: %s" %(self.specie, self.name, self.SNPTypes)
