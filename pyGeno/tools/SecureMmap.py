import mmap

class SecureMmap:
	"""In a normal mmap, modifying the string modifies the file. This is a mmap with write protection"""
	
	def __init__(self, filename, enableWrite = False) :
		
		self.enableWrite = enableWrite
		self.filename = filename
		self.name = filename
		
		f = open(filename, 'r+b')
		self.data = mmap.mmap(f.fileno(), 0)
	
	def forceSet(self, x1, v) :
		"""Forces modification even if the mmap is write protected"""
		self.data[x1] = v
	
	def __getitem__(self, i):
		return self.data[i]
	
	def __setitem__(self, i, v) :
		if self.enableWrite :
			raise IOError("Secure mmap is write protected")
		else :
			self.data[i] = v

	def __str__(self) :
		return "secure mmap, file: %s, writing enabled : %s" % (self.filename, str(self.enableWrite))

	def __len__(self) :
		return len(self.data)
