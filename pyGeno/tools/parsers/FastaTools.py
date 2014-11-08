import os

class FastaFile(object) :
	"""
	Represents a whole Fasta file::
		
		#reading
		f = FastaFile()
		f.parseFile('hop.fasta')
		for line in f :
			print line
		
		#writing
		f = FastaFile()
		f.add(">prot1", "MLPADEV")
		f.save('myFasta.fasta')
	"""
	def __init__(self, fil = None) :
		self.reset()
		if fil != None :
			self.parseFile(fil)
	
	def reset(self) :
		"""Erases everything"""
		self.data = []
		self.currentPos = 0
	
	def parseStr(self, st) :
		"""Parses a string"""
		self.data = st.split('>')[1:]

	def parseFile(self, fil) :
		"""Opens a file and parses it"""
		f = open(fil)
		self.parseStr(f.read())
		f.close()

	def __splitLine(self, li) :
		if len(self.data[li]) != 2 :
			self.data[li] = self.data[li].replace('\r', '\n')
			self.data[li] = self.data[li].replace('\n\n', '\n')
			l = self.data[li].split('\n')
			header = '>'+l[0]
			data = ''.join(l[1:])
			self.data[li] = (header, data)

	def get(self, i) :
		"""returns the ith entry"""
		self.__splitLine(i)
		return self.data[i]
		
	def add(self, header, data) :
		"""appends a new entry to the file"""
		if header[0] != '>' :
			self.data.append(('>'+header, data))
		else :
			self.data.append((header, data))
	
	def save(self, filePath) :
		"""saves the file into filePath"""
		f = open(filePath, 'w')
		f.write(self.make())
		f.close()
	
	def toStr(self) :
		"""returns a string version of self"""
		st = ""
		for d in self.data :
			st += "%s\n%s\n" % (d[0], d[1]) 
	
		return st
		
	def __iter__(self) :
		self.currentPos = 0
		return self
	
	def next(self) :
		#self to call getitem, and split he line if necessary
		return self[self.currentPos]

	def __getitem__(self, i) :
		"""returns the ith entry"""
		return self.get(i)
	
	def __setitem__(self, i, v) :
		"""sets the value of the ith entry"""
		if len(v) != 2: 
			raise TypeError("v must have a len of 2 : (header, data)")
			
		self.data[i] = v
		
	def __len__(self) :
		"""returns the number of entries"""
		return len(self.data)
