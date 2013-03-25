import os

toStr = lambda x:str(x)

def removeDuplicates(fIn, fOut) :
	f = open(fIn)
	legend = f.readline()
	
	data = ''
	h = {}
	h[legend] = 0
	
	lines = f.readlines()
	for l in lines :
		if not h.has_key(l) :
			h[l] = 0
			data += l
			
	f.close()
	f = open(fOut, 'w')
	f.write(legend+data)
	f.close()

def catCSV(folder, output, removeDups = False) :

	strCmd = r"""cat %s/*.csv > %s""" %(folder, output)
	os.system(strCmd)

	if removeDups :
		removeDuplicates(output, output)
	
class CSVFile :
	
	def __init__(self, legend = [], separator = ';') :
		
		self.legend = {}
		if type(legend).__name__ == 'list':
			for i in range(len(legend)) :
				self.legend[legend[i].lower()] = i
		elif type(legend).__name__ == 'dict':
			for k in legend :
				self.legend[k.lower()] = legend[k]
	
		self.header = ''
		self.separator = separator
		self.data = []
		self.notSplitted = []
		self.currentPos = 0
		
	def parse(self, file, separator = ';', stringSeparator = '"') :
		#print "!Champs avec "" mal geres!"
		f = open(file)
		self.data = f.readlines()
		
		f.close()
		self.separator = separator
		self.stringSeparator = stringSeparator
		self.legend = {}
		i = 0
	
		for c in self.data[0].lower().replace(stringSeparator, '').split(separator) :
			legendElement = c.strip()
			if legendElement not in self.legend :
				self.legend[legendElement] = i
			
			i+=1
	
		self.data = self.data[1:]
		self.notSplitted = range(len(self.data))
		
		#self.separtorInStrings = separtorInStrings
		for i in range(len(self.data)) :
			self.__splitLine(i)
		
	def __splitLine(self, l) :
		if  self.data[l] != '' :
			self.data[l] = self.data[l].replace(self.stringSeparator, '')
			self.data[l] = self.data[l].strip().split(self.separator)

			self.notSplitted.remove(l)
	
	def setHeader(self, v) :
		self.header = v
		
	def get(self, line, elt) :
		return self.data[line][self.legend[elt.lower()]].replace(self.stringSeparator, '')
		
	def set(self, line, elt, val) :
		if line in self.notSplitted :
			self.__splitLine(line)

		self.data[line][self.legend[elt.lower()]] = val
		
	def addLine(self, listOrStr = None) :
		if listOrStr != None :
			if type(listOrStr).__name__ == 'str':
				self.notSplitted.append(len(self.data))
				if listOrStr[-1] != '\n' :
					self.data.append(listOrStr + '\n')
				return len(self.data) -1
			elif type(listOrStr).__name__ == 'list':
				l = self.addEmptyLine()
				for i in range(len(listOrStr)) :
					self.data[l][i] = listOrStr[i]
				return len(self.data) -1
		return self.addEmptyLine()
	
	def addEmptyLine(self) :
		l = range(len(self.legend))
		for i in range(len(self.legend)) :
			l[i] =  ''
		self.data.append(l)
		return len(self.data) -1
	
	def addCategoryToLegend(self, category):
		self.legend[category] = len(self.legend)
		for l in self.data :
			l.insert(self.legend[category], '')
			
	def save(self, filePath) :
		f = open(filePath, 'w')
		f.write(self.make())
		f.close()

	def make(self) :
		s = ''
		if len(self.data) > 0 :
			for i in range(len(self.data)) :
				if i in self.notSplitted :
					s += self.data[i]
				else :
					self.data[i] = map(toStr, self.data[i])
					s += self.separator.join(self.data[i]) + '\n'
		
		if self.header != '' :
			return self.header+'\n'+self.separator.join(self.getLegend()) + '\n' + s
		return self.separator.join(self.getLegend()) + '\n' + s
		
	def save(self, filename):
		f = open(filename, 'w')
		f.write(self.make())
		f.close()

	def getLegend(self) :
		sortedLegend = range(len(self.legend))
		for k in self.legend.keys() :
			sortedLegend[self.legend[k]] = k
		return sortedLegend
		
	def __iter__(self) :
		self.currentPos = 0
		return self
	
	def next(self) :
		return self[self.currentPos]
	
	def __getitem__(self, i) :
		if i in self.notSplitted :
			self.__splitLine(i)
		
		return self.data[i]

	def __len__(self) :
		return len(self.data)
