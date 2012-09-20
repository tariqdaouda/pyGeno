class NucFastaConverter:
	"""
	converts to pyGeno format
	ex : N = A/T/C/G
	"""
	
	def __init__(self):#, filePath):
		#f = open(filePath, 'r')
		#self.data = f.read()
		#f.close()
		
		self.format = {'A' : 'A', 'T' : 'T', 'C' : 'C', 'G': 'G',
						'U' : 'T', 'R' : 'A/G', 'Y' : 'C/T', 'M': 'A/C',
						'K' : 'T/G', 'W' : 'A/T', 'S' : 'C/G', 'B':'C/G/T',
						'D' : 'A/G/T', 'H' : 'A/C/T', 'V' : 'A/C/G', 'N': 'A/C/G/T', 
						'X' : 'X', '-' : '-'
				}

		self.format2 = {}
		for k in self.format.keys() :
			self.format2[self.format[k]] = k

	def FastaToPyGeno(self, sequence) :
		s = sequence
		for k in self.format.keys() :
			s = s.replace(k, self.format[k])
		
		return s
	
	def pyGenoToFasta(self, sequence) :
		s = sequence
		for k in self.format2.keys() :
			s = s.replace(k, self.format2[k])
		
		return s