import os, types

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

def catCSVs(folder, output, removeDups = False) :

	strCmd = r"""cat %s/*.csv > %s""" %(folder, output)
	os.system(strCmd)

	if removeDups :
		removeDuplicates(output, output)

def joinCSVs(csvFilePaths, column, output, separator = ';') :
	"csvFilePaths should be an iterable"
	
	res = ''
	legend = []
	csvs = []
	
	for f in csvFilePaths :
		c = CSVFile()
		c.parse(f)
		csvs.append(c)
		legend.append(separator.join(c.legend.keys()))
	legend = separator.join(legend)
	
	lines = []
	for i in range(len(csvs[0])) :
		val = csvs[0].get(i, column)
		line = separator.join(csvs[0][i])
		for c in csvs[1:] :
			for j in range(len(c)) :
				if val == c.get(j, column) :
					line += separator + separator.join(c[j])
					
		lines.append( line )
	print len(lines)
	res = legend + '\n' + '\n'.join(lines)
	
	f = open(output, 'w')
	f.write(res)
	f.close()
	
	return res

class CSVEntry(object) :
	def __init__(self, csvFile, lineNumber = None) :
		
		self.csvFile = csvFile
		self.data = []
		if lineNumber != None :
			self.lineNumber = lineNumber
			tmpL = csvFile.lines[lineNumber].replace('\r', '\n').replace('\n', '')
			tmpData = tmpL.split(csvFile.separator)

			tmpDatum = []
			for d in tmpData :
				sd = d.strip()
				if len(tmpDatum) > 0 or (len(sd) > 0 and sd[0] == csvFile.stringSeparator) :
					tmpDatum.append(sd)
				
					if len(sd) > 0 and sd[-1] == csvFile.stringSeparator :
						self.data.append(csvFile.separator.join(tmpDatum))
						tmpDatum = []
				else :
					self.data.append(sd) 
		else :
			self.lineNumber = len(csvFile)
			for i in range(len(self.csvFile.legend)) :
				self.data.append('')

	def __getitem__(self, key) :
		return self.data[self.csvFile.legend[key.lower()]]

	def __setitem__(self, key, value) :
		self.data[self.csvFile.legend[key.lower()]] = str(value)
	
	def __repr__(self) :
		return "<line %d: %s>" %(self.lineNumber, str(self.data))
		
	def __str__(self) :
		return self.csvFile.separator.join(self.data)
	
class CSVFile(object) :

	def __init__(self, legend = [], separator = ',') :
		"Legend can either be a list of a dict {field : column number}"
		
		self.legend = {}
		if type(legend) is types.ListType :
			for i in range(len(legend)) :
				self.legend[legend[i].lower()] = i
			self.strLegend = separator.join(legend)
			
		elif type(legend) is types.DictType :
			self.strLegend = []
			for k in legend :
				self.legend[k.lower()] = legend[k]
				self.strLegend.insert(legend[k], k.lower())
			self.strLegend = separator.join(self.strLegend)
		
		self.filename = ""
		self.lines = []	
		self.separator = separator
		self.currentPos = -1
	
	def parse(self, filePath, separator = ',', stringSeparator = '"') :
		"Parses a CSV on disc"
		
		self.filename = filePath
		f = open(filePath)
		self.lines = f.readlines()
		f.close()
		
		self.separator = separator
		self.stringSeparator = stringSeparator
		self.legend = {}
		
		i = 0
		for c in self.lines[0].lower().replace(stringSeparator, '').split(separator) :
			legendElement = c.strip()
			if legendElement not in self.legend :
				self.legend[legendElement] = i
			i+=1
	
		self.strLegend = self.lines[0].replace('\r', '\n').replace('\n', '')
		self.lines = self.lines[1:]
		
	def _developLine(self, line) :
		self.lines[line] = CSVEntry(self, line)
	
	def get(self, line, key) :
		if self.lines[line].__class__ is not CSVEntry :
			self._developLine(line)
		
		return self.lines[line][key]

	def set(self, line, key, val) :
		if self.lines[line].__class__ is not CSVEntry :
			self._developLine(line)
		self.lines[line][key] = val
	
	def newLine(self) :
		"Appends an empty line to the end of the CSV and returns it"
		self.lines.append(CSVEntry(self))
		return self.lines[-1]
	
	def insertLine(self, i) :
		"Inserts an empty line at position i and returns it"
		self.data.insert(i, CSVEntry(self))
		return self.lines[i]
	
	def save(self, filePath) :
		"save the csv to a file"
		self.filename = filePath
		f = open(filePath, 'w')
		f.write(self.toStr())
		f.close()

	def toStr(self) :
		s = [self.strLegend]
		for l in self.lines :
			s.append(str(l))
		return '\n'.join(s)
	
	def __iter__(self) :
		self.currentPos = -1
		return self
	
	def next(self) :
		self.currentPos += 1
		if self.currentPos >= len(self) :
			raise StopIteration
		return self[self.currentPos]
	
	def __getitem__(self, line) :
		if self.lines[line].__class__ is not CSVEntry :
			self._developLine(line)
		return self.lines[line]

	def __len__(self) :
		return len(self.lines)
