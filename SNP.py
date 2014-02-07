import configuration as conf

from pyGenoObject import *
import rabaDB.fields as rf

from tools import UsefulFunctions as uf
from exceptions import *

class SNPMaster(Raba) :
	'This object keeps trac'
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	specie = rf.Primitive()
	SNPType = rf.Primitive()
	setName = rf.Primitive()

class SNP_INDEL(pyGenoObject) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	specie = rf.Primitive()
	setName = rf.Primitive()
	chromosomeNumber = rf.Primitive()

	start = rf.Primitive()
	end = rf.Primitive()
	alleles = rf.Primitive()

class CasavaSNP(SNP_INDEL) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	bcalls_used = rf.Primitive()
	bcalls_filt = rf.Primitive()
	ref = rf.Primitive()
	QSNP = rf.Primitive()
	max_gt = rf.Primitive()
	Qmax_gt = rf.Primitive()
	max_gt_poly_site = rf.Primitive()
	Qmax_gt_poly_site = rf.Primitive()
	A_used = rf.Primitive()
	C_used = rf.Primitive()
	G_used = rf.Primitive()
	T_used = rf.Primitive()

	def __init__(self) :
		self.max_gt = self.alleles
		self.pos = self.start

	def _curate(self) :
		self.alleles = self.max_gt
		self.start = self.pos
