import os

toStr = lambda x:str(x)

class GTFFile :
	
	def __init__(self) :
		
		self.data = []
		self.splittedFlag = []
		self.currentPos = 0
		#self.openedLines = []
		leg = ['chromosome', 'coding_type', 'region_type', 'x1', 'x2', '.1', 'strand', '.2',
						'gene_id', 'transcript_id','exon_number', 'gene_name', 'transcript_name',
						'protein_id']
		self.legend ={}
		i = 0
		for l in leg:
			self.legend[l] = i
			i+=1
		
	def parseStr(self, str) :
		self.data = str.split('\n')
		if self.data[-1] == '':
			del self.data[-1]
		
		self.splittedFlag = range(len(self.data))

	def parseFile(self, file) :
		f = open(file)
		self.data = f.readlines()
		f.close()		
		self.splittedFlag = range(len(self.data))
		
	def __splitLine(self, l) :
		if self.splittedFlag[l] != 'split' :
			
			self.data[l] = self.data[l].replace(' ', '').replace('\t', ';').replace('"', '')
			for v in ['gene_id','transcript_id','exon_number','gene_name','transcript_name','protein_id']:
				self.data[l] = self.data[l].replace(v, '')
				
			self.data[l] = self.data[l].split(';')
			if self.data[l][-1] == '':
				del self.data[l][-1]
			
			self.splittedFlag[l] = 'split'
		
	def getElement(self, line, elt) :
		self.__splitLine(line)
		try :
			return self.data[line][self.legend[elt.lower()]]
		except IndexError :
			return ''
			
	def getField(self, line, elt) :
		return self. getElement(line, elt)
	
	def get(self, line, elt) :
		return self. getElement(line, elt)
		
	def setElement(self, line, elt, val) :
		self.__splitLine(line)
		self.data[line][self.legend[elt.lower()]] = val
	
	def setField(self, line, elt, val) :
		self. setElement(line, elt, val)
	
	def set(self, line, elt, val) :
		self. setElement(line, elt, val)
		
	def setElementsToLine(self, line, elts) :
		self.__splitLine(line)
		self.data[line].extend(elts)
	
	def addLine(self, listy) :
		self.data.append(listy)
		return len(self.data) -1
	
	"""	
	def save(self, filePath) :
		f = open(filePath, 'w')
		f.write(self.make())
		f.close()
	
	def toStr(self) :
		return self.make()
	
	def build(self) :
		return self.make()
		
	def make(self) :
		#print self.legend
		s = ''
		if len(self.data) > 0 :
			for i in range(len(self.data)) :
				if i in self.notSplitted :
					s += self.data[i]
				else :
					self.data[i] = map(toStr, self.data[i])
					s += self.separator.join(self.data[i]) + '\n'
		
			
		return self.separator.join(self.getLegend()) + '\n' + s
	
	def getLegend(self) :
		sortedLegend = range(len(self.legend))
		for k in self.legend.keys() :
			sortedLegend[self.legend[k]] = k
		return sortedLegend
		
	def __iter__(self) :
		self.currentPos = 0
		return self
	
	def next(self) :
		#self to call getitem, and split he line if necessary
		return self[self.currentPos]
	"""
	def __getitem__(self, i) :
		#print self.data[i], i in self.notSplitted
		if i in self.notSplitted :
			self.__splitLine(i)
		ret = {}
		for k in self.legend:
			ret[k] = self.data[i][self.legend[k]]
		
		return ret

	def __len__(self) :
		return len(self.data)
