import gzip

class GTFFile :
	"""This is a simple GTF2.2 (Revised Ensembl GTF) parser, see http://mblab.wustl.edu/GTF22.html for more infos"""
	def __init__(self) :
		
		self.data = []
		self.splittedFlag = [] # < 0 means line has already been processed
		self.currentPos = 0
		self.legend = {'seqname' : 0, 'source' : 1, 'feature' : 2, 'start' : 3, 'end' : 4, 'score' : 5, 'strand' : 6, 'frame' : 7, 'attributes' : 8}
	
	def parseStr(self, string) :
		self.data = string.split('\n')
		if self.data[-1] == '':
			del self.data[-1]
		
		self.splittedFlag = range(len(self.data))

	def parseFile(self, fil, gziped = False) :
		if gziped : 
			f = gzip.open(fil)
		else :
			f = open(fil)
		
		self.data = f.readlines()
		f.close()
		if self.data[-1] == '':
			del self.data[-1]
		
		self.splittedFlag = range(len(self.data))
		
	def __splitLine(self, l) :
		#print '---<', l, self.data[l]
		if self.data[l][0] == '#' :
			return False
		
		if self.splittedFlag[l] >= 0 :
			self.data[l] = self.data[l][:-2].split('\t') #-2 remove ';\n'
			proto_atts = self.data[l][self.legend['attributes']].strip().split('; ')
			atts = {}
			for a in proto_atts :
				sa = a.split(' ')
				atts[sa[0]] = sa[1].replace('"', '')	
			self.data[l][self.legend['attributes']] = atts
			
			self.splittedFlag[l] = -l-1
		return True
	
	def get(self, line, elmt) :
		i = line
		while self.__splitLine(i) == None :
			i += 1
		
		try :
			elmtId = self.legend[elmt]
			#print 'iyoiuy', line, elmt, self.data[i]
			return self.data[i][elmtId]
		except KeyError :
			try :
				return self.data[i][self.legend['attributes']][elmt]
			except KeyError :
				raise KeyError("Line %d does not have an element %s. Line : %s" %(i, elmt, self.data[i]))
		
	def __getitem__(self, i) :
		if i in self.notSplitted :
			self.__splitLine(i)
		ret = {}
		for k in self.legend:
			ret[k] = self.data[i][self.legend[k]]
		
		return ret

	def __len__(self) :
		return len(self.data)
