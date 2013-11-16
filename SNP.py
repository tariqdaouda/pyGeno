import configuration as conf

from rabaDB.setup import *
RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
from rabaDB.Raba import *
import rabaDB.fields as rf

from tools import UsefulFunctions as uf
from exceptions import *

class Casava_SNP(Raba) :
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	pos = rf.PrimitiveField()
	bcalls_used = rf.PrimitiveField()
	bcalls_filt = rf.PrimitiveField()
	ref = rf.PrimitiveField()
	QSNP = rf.PrimitiveField()
	max_gt = rf.PrimitiveField()
	max_gt_poly_site = rf.PrimitiveField()
	Qmax_gt_poly_site = rf.PrimitiveField()
	A_used = rf.PrimitiveField()
	C_used = rf.PrimitiveField()
	G_used = rf.PrimitiveField()
	T_used = rf.PrimitiveField()
	
	genome = rf.RabaObjectField('Genome')
	chromosome = rf.RabaObjectField('Chromosome')
	
	def __init__(self, *args, **fieldsSet) :
		Raba.__init__(self, **fieldsSet) 

class dbSNP_SNP(Raba) :
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	version = rf.PrimitiveField()
	
	rsid = rf.PrimitiveField()
	pos = rf.PrimitiveField()
	type = rf.PrimitiveField()
	alleles = rf.PrimitiveField()
	validated = rf.PrimitiveField()
	
	assembly = rf.PrimitiveField()
	original_orientation = rf.PrimitiveField()
	maf_allele = rf.PrimitiveField()
	maf_count = rf.PrimitiveField()
	maf = rf.PrimitiveField()
	het = rf.PrimitiveField()
	se_het = rf.PrimitiveField()
	
	loc = rf.RabaListField()
	
	specie = rf.PrimitiveField()
	chromosomeNumber = rf.PrimitiveField()
	
	def __init__(self, *args, **fieldsSet) :
		Raba.__init__(self, **fieldsSet) 

class dbSNP_SNPLOC(Raba) :
	
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	allele = rf.PrimitiveField()
	fxn_class = rf.PrimitiveField()
	gene = rf.PrimitiveField()
	residue = rf.PrimitiveField()
	
	snp = rf.RabaObjectField('dbSNP_SNP')
	
	def __init__(self, *args, **fieldsSet) :
		Raba.__init__(self, **fieldsSet)
		
"""
class CasavaSNP(SNP) :
	def __init__(self, line) :
		SNP.__init__(self, line)
		self.formatDecriptionFile = 'casavaSNP_FormatDescription.txt'
		self.__make(line.split(';'))
		
	def __make(self, sl) :
		if len(sl) != 13 :
			raise SNPError("a Casava SNP should have 13 fields, got: %d" % len(sl), sl)
		
		try :
			self.values = {}
			self.values['pos'] = int(sl[0])
			self.values['bcalls_used'] = sl[1]
			self.values['bcalls_filt'] = sl[2]
			self.values['ref'] = sl[3]
			self.values['QSNP'] = int(sl[4])
			self.values['max_gt'] = uf.getPolymorphicNucleotide(sl[5])
			self.values['Qmax_gt'] = int(sl[6])
			self.values['max_gt_poly_site'] = sl[7]
			self.values['Qmax_gt_poly_site'] = int(sl[8])
			self.values['A_used'] = int(sl[9])
			self.values['C_used'] = int(sl[10])
			self.values['G_used'] = int(sl[11])
			self.values['T_used'] = int(sl[12])
		except :
			raise  SNPError("Unkown problem, here's the line: %s" % sl)

class dbSNP(SNP) :
	def __init__(self, line) :
		SNP.__init__(self, line)
		self.formatDecriptionFile = 'dbSNP_FormatDescription.txt'
		self.__make(line.split(';'))
	
	def __make(self, sl) :
		#print '---<', sl
		if len(sl) != 14 :
			raise SNPError("a dbSNP SNP should have 14 fields, got: %d" % len(sl), sl)
		
		try :
			self.values = {}
			self.values['pos'] = int(sl[0])
			self.values['chr'] = sl[1]
			self.values['rs'] = int(sl[2])
			self.values['type'] = sl[3]

			try :
				self.values['alleles'] = uf.getPolymorphicNucleotide(sl[4])
			except uf.UnknownNucleotide :	
				self.values['alleles'] = sl[4]

			self.values['validated'] = (sl[5].upper() == 'YES')
			self.values['assembly'] = sl[6]
			self.values['original_strand'] = sl[7]
			self.values['maf_allele'] = sl[8]
			self.values['maf_count'] = int(float(sl[9]))
			self.values['maf'] = float(sl[10])
			self.values['het'] = float(sl[11])
			self.values['se(het)'] = float(sl[12])
			#print pickle.loads(zlib.decompress(sl[13]))
			self.values['loc'] = pickle.loads(sl[13].replace('/rje3/', '\n').replace('/qte3/', ';'))
		except :
			raise SNPError("Unkown problem, here's the line: %s" % sl)
"""
