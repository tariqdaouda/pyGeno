import os

class FastqEntry(object) :
	"""A single entry in the FastqEntry file"""
	
	def __init__(self, ident = "", seq = "", plus = "", qual = "") :
		self.values = {}
		self.values['identifier'] = ident
		self.values['sequence'] = seq
		self.values['+'] = plus
		self.values['qualities'] = qual
	
	def __getitem__(self, i):
		return self.values[i]
	
	def __setitem__(self, i, v) :
		self.values[i] = v
		
	def __str__(self):
		return "%s\n%s\n%s\n%s" %(self.values['identifier'], self.values['sequence'], self.values['+'], self.values['qualities'])
	
class FastqFile(object) :
	"""
	Represents a whole CSV file::
		
		#reading
		f = FastqFile()
		f.parse('hop.csv')
		for line in f :
			print line['sequence']
		
		#writing, legend can either be a list of a dict {field : column number}
		f = CSVFile(legend = ['name', 'email'])
		l = f.newLine()
		l['name'] = 'toto'
		l['email'] = "hop@gmail.com"
		f.save('myCSV.csv')
		
	"""
	
	def __init__(self, fil = None) :
		self.reset()
		if fil != None :
			self.parseFile(fil)
	
	def reset(self) :
		"""Frees the file"""
		self.data = []
		self.currentPos = 0
	
	def parseStr(self, st) :
		"""Parses a string"""
		self.data = st.replace('\r', '\n')
		self.data = self.data.replace('\n\n', '\n')
		self.data = self.data.split('\n')

	def parseFile(self, fil) :
		"""Parses a file on disc"""
		f = open(fil)
		self.parseStr(f.read())
		f.close()		
		
	def __splitEntry(self, li) :
		try :
			self.data['+'] 
			return self.data
		except :
			self.data[li] = FastqEntry(self.data[li], self.data[li+1], self.data[li+2], self.data[li+3])
			
	def get(self, li) :
		"""returns the ith entry"""
		i = li*4
		self.__splitEntry(i)
		return self.data[i]
	
	def newEntry(self, ident = "", seq = "", plus = "", qual = "") :
		"""Appends an empty entry at the end of the CSV and returns it"""
		e = FastqEntry()
		self.data.append(e)
		return e
	
	def add(self, fastqEntry) :
		"""appends an entry to self"""
		self.data.append(fastqEntry)
		
	def save(self, filePath) :
		f = open(filePath, 'w')
		f.write(self.make())
		f.close()
			
	def toStr(self) :
		st = ""
		for d in self.data :
			st += "%s\n%s" % (d[0], d[1]) 
	
		return st
		
	def __iter__(self) :
		self.currentPos = 0
		return self
	
	def next(self) :
		#self to call getitem, and split he line if necessary
		i = self.currentPos +1
		#print i-1, self.currentPos
		if i > len(self) :
			raise StopIteration()
			
		self.currentPos = i
		return self[self.currentPos-1]

	def __getitem__(self, i) :
		"""returns the ith entry"""
		return self.get(i)
	
	def __setitem__(self, i, v) :
		"""sets the ith entry"""
		if len(v) != 2: 
			raise TypeError("v must have a len of 2 : (header, data)")
			
		self.data[i] = v
		
	def __len__(self) :
		return len(self.data)/4
