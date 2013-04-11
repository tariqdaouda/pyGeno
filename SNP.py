import configuration as conf
from tools import UsefulFunctions as uf

class SNPError(Exception) :
	def __init__(self, message, line):
		self.message = message
		self.line = line
	
	def __str__(self):
		return "Invalid SNP, %s" % self.message

class SNPFile :
	"""This represent a file containing a list of snps sorted by position, the position must be in the first column"""
	def __init__(self, filePath, SNPObj) :
		self.filePath = filePath
		f = open(filePath)
		self.lines = f.readlines()
		f.close()
		if self.lines[0][:2] == '//' :
			j = 0
			for i in range(len(self.lines)) :
				if self.lines[i][:6] =='//pos;' :
					break
				j += 1
				
			self.lines = self.lines[j+1:]
			
		if len(self.lines) < 1 :
			raise ValueError("file %s is empty" % filePath)

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

	def findSnpsInRange(self, x1, x2 = None) :
		"""X1 inclusive, x2 exc"""
		if x2 == None or x1 == x2:
			xx1 = x1
			xx2 = x1+1
		else :
			if x1 < x2 :
				xx1, xx2 = x1, x2
			else :
				xx1, xx2 = x2, x1
	
		l1, val1 = self.__findSnp(xx1)
		l2, val2 = self.__findSnp(xx2)
		
		ret = []
		for l in range(l1, l2+1) :
			snp = self.__getSNP(l)
			if  xx1<= snp['pos'] and snp['pos'] < xx2 :
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
	
	def __repr__(self):
		return str(self)
		
	def __eq__(self, s) :
		return self.values == s.values
	
	def printFormatDescription() :
		f = open(conf.DATA_PATH+'/'+self.formatDecriptionFile)
		s = f.read()
		f.close()
		print s
	
	def __hash__(self):
		return hash(self.__class__.__name__) ^ hash(self.values)
	
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
			self.values['max_gt-poly_site'] = sl[7]
			self.values['Qmax_gt-poly_site'] = int(sl[8])
			self.values['A_used'] = int(sl[9])
			self.values['C_used'] = int(sl[10])
			self.values['G_used'] = int(sl[11])
			self.values['T_used'] = int(sl[12])
		except :
			raise SNPError("Unkown problem, see self.line")

class dbSNP(SNP) :
	def __init__(self, line) :
		SNP.__init__(self, line)
		self.formatDecriptionFile = 'dbSNP_FormatDescription.txt'
		self.__make(line.split(';'))
	
	def __make(self, sl) :
		#print sl
		if len(sl) != 13 :
			raise SNPError("a dbSNP SNP should have 13 fields, got: %d" % len(sl), sl)
		
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
		except :
			raise SNPError("Unkown problem, see self.line")

