import configuration as conf

from pyGenoObject import *
import rabaDB.fields as rf

from tools import UsefulFunctions as uf
from exceptions import *

class CasavaSNP(pyGenoObject) :

	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	pos = rf.Primitive()
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

	genome = rf.RabaObject('Genome')
	chromosome = rf.RabaObject('Chromosome')

	def __init__(self) :
		pass

	def _curate(self) :
		pass
