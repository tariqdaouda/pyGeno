import configuration as conf

from rabaDB.setup import *
RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
from rabaDB.Raba import *
import rabaDB.fields as rf

from tools import UsefulFunctions as uf
from exceptions import *

class Casava_SNP(Raba) :
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	pos = rf.Primitive()
	bcalls_used = rf.Primitive()
	bcalls_filt = rf.Primitive()
	ref = rf.Primitive()
	QSNP = rf.Primitive()
	max_gt = rf.Primitive()
	max_gt_poly_site = rf.Primitive()
	Qmax_gt_poly_site = rf.Primitive()
	A_used = rf.Primitive()
	C_used = rf.Primitive()
	G_used = rf.Primitive()
	T_used = rf.Primitive()
	
	genome = rf.RabaObject('Genome')
	chromosome = rf.RabaObject('Chromosome')
	
	def __init__(self, *args, **fieldsSet) :
		pass

class dbSNP_SNP(Raba) :
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	version = rf.Primitive()
	
	rsId = rf.Primitive()
	pos = rf.Primitive()
	type = rf.Primitive()
	alleles = rf.Primitive()
	validated = rf.Primitive()
	
	assembly = rf.Primitive()
	original_orientation = rf.Primitive()
	maf_allele = rf.Primitive()
	maf_count = rf.Primitive()
	maf = rf.Primitive()
	het = rf.Primitive()
	se_het = rf.Primitive()
	
	loc = rf.RList()
	
	specie = rf.Primitive()
	chromosomeNumber = rf.Primitive()
	
	def __init__(self, *args, **fieldsSet) :
		pass

class dbSNP_SNPLOC(Raba) :
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	allele = rf.Primitive()
	fxn_class = rf.Primitive()
	gene = rf.Primitive()
	residue = rf.Primitive()
	
	snp = rf.RabaObject('dbSNP_SNP')
	
	def __init__(self, *args, **fieldsSet) :
		pass
