import configuration as conf
from tools import UsefulFunctions as uf

class SNPFile :
	"""This represent a file containing a list of snps sorted by position, the position must be in the first column"""
	def __init__(self, filePath, SNPObj) :
		f = open(filePath)
		self.lines = f.readlines()
		f.close()
		self.snps = {}
		self.SNPObj = SNPObj
	def __findSnp(self, x1):
		r1 = 0
		r2 = len(self)-1
		while (r1 <= r2) :
			pos = (r1+r2)/2
			sl = self.lines[pos].split(';')
			val = int(sl[0])

			if val == x1 :
				return (pos, val)
			elif x1 < val :
				r2 = pos -1
			else :
				r1 = pos +1

		return (pos, val)

	def findSnpsInRange(self, x1, x2) :
		"""X1 inclusive, x2 exc"""
		l1, val1 = self.__findSnp(x1)
		l2, val2 = self.__findSnp(x2)
		
		ret = []
		for l in range(l1, l2+1) :
			snp = self.__getSNP(l)
			if  x1<= snp['pos'] and snp['pos'] <x2 :
				ret.append(snp)

		return ret
	
	def findSnp(self, x1) :
		l1, val1 = self.__findSnp(x1)

		if val1 == x1 :
			return self.__getSNP(l1)
		else :
			return None
			
	def __getSNP(self, lineNumber) :
		if lineNumber > len(self) or lineNumber < 0:
			return None
		
		try :
			return self.snps[lineNumber]
		except KeyError :
			#self.snps[lineNumber] = CasavaSNP(self.lines[lineNumber])
			self.snps[lineNumber] = self.SNPObj(self.lines[lineNumber])
			return self.snps[lineNumber]
	
	def __len__(self) :
		return len(self.lines)

class SNP :
	def __init__(self, line) :
		self.reset(line)
	
	def reset(self, line) :
		"""An abstract class representing a SNP, see __make()"""
		self.values = {}
		self.formatDecriptionFile = ''
		
	def __make(self, sl) :
		"""This function should be rewritten in child, i fills the self.value dictionary from a list"""
		raise TypeError("SNP Class is abstact and should never get instanciated")
		
	def __getitem__(self, i) :
		return self.values[i]
	
	def __setitem__(self, i, v) :
		self.values[i] = v
		
	def __str__(self) :
		return str(self.values)
	
	def __eq__(self, s) :
		return self.values == s.values
	
	def printFormatDescription() :
		f = open(conf.DATA_PATH+'/'+self.formatDecriptionFile)
		s = f.read()
		f.close()
		print s

class CasavaSNP(SNP) :
	def __init__(self, line) :
		SNP.__init__(self, line)
		self.formatDecriptionFile = 'casavaSNP_FormatDescription.txt'
		self.__make(line.split(';'))
		
	def __make(self, sl) :		
		self.values = {}
		self.values['pos'] = int(sl[0])
		self.values['bcalls_used'] = sl[1]
		self.values['bcalls_filt'] = sl[2]
		self.values['ref'] = sl[3]
		self.values['QSNP'] = int(sl[4])
		self.values['max_gt'] = uf.getPolymorphicNucleotide(sl[5])
		self.values['Qmax_gt'] = int(sl[6])
		self.values['max_gt-poly_site'] = sl[7]
		self.values['Qmax_gt-poly_site'] = int(sl[8])
		self.values['A_used'] = int(sl[9])
		self.values['C_used'] = int(sl[10])
		self.values['G_used'] = int(sl[11])
		self.values['T_used'] = int(sl[12])


class dbSNP(SNP) :
	def __init__(self, line) :
		SNP.__init__(self, line)
		self.formatDecriptionFile = 'dbSNP_FormatDescription.txt'
		self.__make(line.split(';'))
		
	def __make(self, sl) :
		self.values = {}
		self.values['pos'] = int(sl[0])
		self.values['mapwgt'] = int(sl[1])
		self.values['snp_type'] = int(sl[2])
		self.values['chr_hits'] = int(sl[3])
		self.values['ctg_hits'] = int(sl[4])
		self.values['total_hits'] = int(sl[5])
		self.values['chr'] = sl[6]
		self.values['ctg_acc'] = sl[7]
		self.values['ctg_ver'] = int(sl[8])
		self.values['ctg_ID'] = sl[9]
		self.values['ctg_pos'] = int(sl[10])
		self.values['rs'] = int(sl[11])
		self.values['local_loci'] = int(sl[12])
		self.values['avg_het'] = int(sl[13])
		self.values['s.e._het'] = int(sl[14])
		self.values['max_prob'] = int(sl[15])
		self.values['validated'] = int(sl[16])
		self.values['genotypes'] = int(sl[17])
		self.values['links_out'] = int(sl[18])
		self.values['orig_build'] = sl[19]
		self.values['upd_build'] = int(sl[20])

