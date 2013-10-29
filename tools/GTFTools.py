class GTFFile :
	"""This is a simple GTF2.2 (Revised Ensembl GTF) parser, see http://mblab.wustl.edu/GTF22.html for more infos"""
	def __init__(self) :
		
		self.data = []
		self.splittedFlag = [] # < 0 means line has already been processed
		self.currentPos = 0
		legend = {'seqname' : 0, 'source' : 1, 'feature' : 2, 'start' : 3, 'end' : 4, 'score' : 5, 'frame' : 6, 'attributes' : 7}
	
	def parseStr(self, str) :
		self.data = str.split('\n')
		if self.data[-1] == '':
			del self.data[-1]
		
		self.splittedFlag = range(len(self.data))

	def parseFile(self, file) :
		f = open(file)
		self.data = f.readlines()
		f.close()
		if self.data[-1] == '':
			del self.data[-1]
		
		self.splittedFlag = range(len(self.data))
		
	def __splitLine(self, l) :
		if self.splittedFlag[l] >= 0 :
			self.data[l] = self.data[l].split('\t')
			
			proto_atts = self.data[l][self.legend['attributes']].split('; ')
			atts = {}
			for a in atts :
				sa = a.split(' ')
				atts[sa[0]] = sa[1].replace('"', '')	
			self.data[l][self.legend['attributes']] = atts
			
			self.splittedFlag[l] = -l
		
	def get(self, line, elt) :
		self.__splitLine(line)
		try :
			return self.data[line][elmt]
		except IndexError :
			try :
				att = self.data[line][self.legend['attributes']][elmt]
			except IndexError :
				None
		
	def __getitem__(self, i) :
		if i in self.notSplitted :
			self.__splitLine(i)
		ret = {}
		for k in self.legend:
			ret[k] = self.data[i][self.legend[k]]
		
		return ret

	def __len__(self) :
		return len(self.data)
