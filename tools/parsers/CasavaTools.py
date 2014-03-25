import UsefulFunctions as uf

class SNPsTxtEntry :

	def __init__(self, line) :
		self.values = {}
		sl = line.replace('\t\t', '\t').split('\t')
		
		self.values['chromosomeNumber'] = sl[0].upper().replace('CHR', '')
		#first column: chro, second first of range (identical to second column)
		self.values['start'] = int(sl[2])
		self.values['end'] = int(sl[2])+1
		self.values['bcalls_used'] = sl[3]
		self.values['bcalls_filt'] = sl[4]
		self.values['ref'] = sl[5]
		self.values['QSNP'] = int(sl[6])
		self.values['alleles'] = uf.getPolymorphicNucleotide(sl[7]) #max_gt
		self.values['Qmax_gt'] = int(sl[8])
		self.values['max_gt_poly_site'] = sl[9]
		self.values['Qmax_gt_poly_site'] = int(sl[10])
		self.values['A_used'] = int(sl[11])
		self.values['C_used'] = int(sl[12])
		self.values['G_used'] = int(sl[13])
		self.values['T_used'] = int(sl[14])
	
	def __getitem__(self, i):
		return self.values[i]
	
	def __setitem__(self, i, v) :
		self.values[i] = v
		
	def __str__(self):
		return str(self.values)
	
class SNPsTxtFile :
	
	def __init__(self, fil = None) :
		self.reset()
		if fil != None :
			self.parseFile(fil)
	
	def reset(self) :
		self.data = []
		self.currentPos = 0
	
	def parseStr(self, st) :
		lines = st.replace('\r', '\n')
		lines = self.data.replace('\n\n', '\n')
		lines = self.data.split('\n')
		
		for l in lines :
			if l[0] != '#'
			self.data.append(SNPsTxtEntry(l))
	
	def parseFile(self, fil) :
		f = open(fil)
		self.parseStr(f.read())
		f.close()		
		
	def get(self, li) :
		return self.data[i]
		
	def save(self, filePath) :
		f = open(filePath, 'w')
		f.write(self.make())
		f.close()
	
	def toStr(self) :
		return self.make()
	
	def build(self) :
		return self.make()
		
	def make(self) :
		st = ""
		for d in self.data :
			st += "%s\n%s" % (d[0], d[1]) 
	
		return st
		
	def __iter__(self) :
		self.currentPos = 0
		return self
	
	def next(self) :
		if self.currentPos >= len(self) :
			raise StopIteration()
		v = self[self.currentPos]
		self.currentPos += 1
		return v

	def __getitem__(self, i) :
		return self.get(i)
	
	def __setitem__(self, i, v) :
		if len(v) != 2: 
			raise TypeError("v must have a len of 2 : (header, data)")
			
		self.data[i] = v
		
	def __len__(self) :
		return len(self.data)
