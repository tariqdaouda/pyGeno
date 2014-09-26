import os

class FastqEntry(object) :
	def __init__(self, ident, seq, plus, qual) :
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
	
	def __init__(self, fil = None) :
		self.reset()
		if fil != None :
			self.parseFile(fil)
	
	def reset(self) :
		self.data = []
		self.currentPos = 0
	
	def parseStr(self, st) :
		self.data = st.replace('\r', '\n')
		self.data = self.data.replace('\n\n', '\n')
		self.data = self.data.split('\n')

	def parseFile(self, fil) :
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
		i = li*4
		self.__splitEntry(i)
		return self.data[i]
		
	def add(self, fastqEntry) :
		self.data.append(fastqEntry)
		
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
		#self to call getitem, and split he line if necessary
		i = self.currentPos +1
		#print i-1, self.currentPos
		if i > len(self) :
			raise StopIteration()
			
		self.currentPos = i
		return self[self.currentPos-1]

	def __getitem__(self, i) :
		return self.get(i)
	
	def __setitem__(self, i, v) :
		if len(v) != 2: 
			raise TypeError("v must have a len of 2 : (header, data)")
			
		self.data[i] = v
		
	def __len__(self) :
		return len(self.data)/4
