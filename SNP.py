import configuration as conf

from pyGenoObjectBases import *
import rabaDB.fields as rf

from tools import UsefulFunctions as uf
from exceptions import *

class SNPMaster(Raba) :
	'This object keeps track of set types'
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	specie = rf.Primitive()
	SNPType = rf.Primitive()
	setName = rf.Primitive()
	_raba_uniques = [('setName',)]

	def _curate(self) :
		self.specie = self.specie.lower()
		self.setName = self.setName.lower()

class SNP_INDEL(pyGenoRabaObject) :
	"All SNPs should inherit from me"
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	_raba_abstract = True # not saved in db

	specie = rf.Primitive()
	setName = rf.Primitive()
	chromosomeNumber = rf.Primitive()

	start = rf.Primitive(default = 3)
	end = rf.Primitive()
	alleles = rf.Primitive()

	def _curate(self) :
		self.specie = self.specie.lower()
		self.setName = self.setName.lower()

class CasavaSNP(SNP_INDEL) :
	"A SNP of Casava"
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE

	bcalls_used = rf.Primitive()
	bcalls_filt = rf.Primitive()
	ref = rf.Primitive()
	QSNP = rf.Primitive()
	Qmax_gt = rf.Primitive()
	max_gt_poly_site = rf.Primitive()
	Qmax_gt_poly_site = rf.Primitive()
	A_used = rf.Primitive()
	C_used = rf.Primitive()
	G_used = rf.Primitive()
	T_used = rf.Primitive()

	def _curate(self) :
		SNP_INDEL._curate(self)

class dbSNPSNP(SNP_INDEL) :
	pass

class TopHatSNP(SNP_INDEL) :
	pass
