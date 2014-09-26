import os

class FastaFile(object) :
	def __init__(self, fil = None) :
		self.reset()
		if fil != None :
			self.parseFile(fil)
	
	def reset(self) :
		self.data = []
		self.currentPos = 0
	
	def parseStr(self, st) :
		self.data = st.split('>')[1:]

	def parseFile(self, file) :
		f = open(file)
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

	def get(self, l) :
		self.__splitLine(l)
		return self.data[l]
		
	def add(self, header, data) :
		if header[0] != '>' :
			self.data.append(('>'+header, data))
		else :
			self.data.append((header, data))
	
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
			st += "%s\n%s\n" % (d[0], d[1]) 
	
		return st
		
	def __iter__(self) :
		self.currentPos = 0
		return self
	
	def next(self) :
		#self to call getitem, and split he line if necessary
		return self[self.currentPos]

	def __getitem__(self, i) :
		return self.get(i)
	
	def __setitem__(self, i, v) :
		if len(v) != 2: 
			raise TypeError("v must have a len of 2 : (header, data)")
			
		self.data[i] = v
		
	def __len__(self) :
		return len(self.data)
