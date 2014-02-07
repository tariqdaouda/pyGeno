import os
import configuration as conf
from pyGenoObject import *

#from rabaDB.setup import *
#RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
#from rabaDB.Raba import *
#from rabaDB.filters import RabaQuery
import rabaDB.fields as rf

class SequenceSNP :
	"This SNP is intended to be embeded in a sequence and has no persistency. It's alleles can be a mix of alleles from SNPs of different sets"
	def __init__(self, alleles, positions) :
		self.alleles = alleles
		self.position = position
		self.sourceSNPs = []

	def addSourceSNP(self, sourceSNP) :
		"persistent SNPs whose alleles have been mixed to make self.alleles"
		self.sourceSNPs.append(sourceSNP)

def defaultSNPsFilter(*SNPs) :
	retSNP = SequenceSNP(None, position = SNPs[0].position)
	alleles = []

	for SNP in SNPs :
		alleles.append(SNP.alleles)
		retSNP.addSNP(SNP)

	retSNP.alleles = ''.join(SNP.alleles)

	return retSNP

class Genome(pyGenoObject) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	name = rf.Primitive()
	specie = rf.Primitive()
	#chromosomes = rf.Relation('Chromosome')
	#dataType = rf.Primitive()

	source = rf.Primitive()
	packageInfos = rf.Primitive()

	#_raba_uniques = [('name', 'specie')]

	def __init__(self, SNPs = None, SNPsFilter = defaultSNPsFilter) :
		self.SNPsSets = SNPs
		self.SNPsFilter = SNPsFilter
		self.SNPTypes = {}

		if SNPs != None :
			f = RabaQuery(SNPMaster, namespace = self._raba_namespace)
			for se in self.SNPsSets :
				f.addFilter(setName = self.SNPsSet, specie = self.specie)

			res = f.run()
			if res == None :
				raise ValueError("There's no set of SNPs that goes by the name of %s for specie %s" % (setName, self.specie))

			for s in res :
				self.SNPTypes[s.setName] = s.SNPType

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

if False :
	m = SequencedGenome(specie = 'human', name = 'ARN_M')
	ref = Genome(specie = 'human', name = 'ARN_M')

	m = Genome(specie = 'human', name = 'GRCH37.73', SNPs = ('ARN_M', 'ARN_R'), snpFilter = fct)
	m = Genome(specie = 'human', SNPs = 'ARN_M')

	m.get(Transcript)[0].sequence
	gene = m.get(Gene)[0]
	print gene.Get(Exon)[0].sequence
