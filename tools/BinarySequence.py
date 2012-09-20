import array, copy
import UsefulFunctions as uf

class BinarySequence :
	
	def __init__(self, sequence, arrayForma, charToBinDict) :
		
		self.polymorphisms = []
		self.defaultSequence = ''
		
		self.forma = arrayForma
		self.charToBin = charToBinDict
		self.sequence = sequence
		
		self.binSequence = self.encode(sequence)
		self.itemsize = self.binSequence.itemsize
		self.typecode = self.binSequence.typecode
		#print 'bin', len(self.sequence), len(self.binSequence)
		
	def encode(self, sequence):
		res = array.array(self.forma.typecode)
		b = 0
		i = 0
		trueI = 0 #not inc in case if poly
		poly = []
		while i < len(sequence)-1:
			b = b | self.forma[self.charToBin[sequence[i]]]
			if sequence[i+1] == '/' :
				if sequence[i] not in poly :
					poly.append(sequence[i])
				i += 2
			else :
				res.append(b)
				b = 0
				if len(poly) > 0 :
					if sequence[i] not in poly :
						poly.append(sequence[i])
					self.polymorphisms.append((trueI, poly))
					poly = []
				self.defaultSequence += sequence[i]
				i += 1
				trueI += 1
				
		if i < len(sequence) :
			b = b | self.forma[self.charToBin[sequence[i]]]
			res.append(b)
			if len(poly) > 0 :
				if sequence[i] not in poly :
					poly.append(sequence[i])
				self.polymorphisms.append((trueI, poly))
				#poly = []
			self.defaultSequence += sequence[i]
		return res
		
	def __testFind(self, arr) :
		if len(arr)  == 0:
			raise TypeError ('binary find, needle is empty')
		if arr.itemsize != self.itemsize :
			raise TypeError ('binary find, both arrays must have same item size, arr: %d, self: %d' % (arr.itemsize, self.itemsize))
	
	def __testBinary(self, arr) :
		if len(arr) != len(self) :
			raise TypeError ('bin operator, both arrays must be of same length')
		if arr.itemsize != self.itemsize :
			raise TypeError ('bin operator, both arrays must have same item size, arr: %d, self: %d' % (arr.itemsize, self.itemsize))
	
	def findPolymorphisms(self, strSeq, strict = False):
		"""If strict = False, it ignores the cases of matching heterozygocity (ex: for a given position i, strSeq[i] = A and self.sequence = 'A/G')
		, if strict = True returns all positions where strSeq differs self's sequence"""
		arr = self.encode(strSeq)
		res = []
		if not strict :
			for i in range(len(arr)+len(self)) :
				if i >= len(arr) or i > len(self) :
					break
				if arr[i] & self[i] == 0:
					res.append(i)
		else :
			for i in range(len(arr)+len(self)) :
				if i >= len(arr) or i > len(self) :
					break
				if arr[i] != self[i] :
					res.append(i)
		return res
		
	def getPolymorphisms(self):
		"returns all polymorphsims in the form of a dict pos => alleles"
		return self.polymorphisms
	
	def getDefaultSequence(self) :
		"returns a version str sequence where only the last allele of each polymorphisms is shown" 
		return self.defaultSequence
	
	def __getSequenceVariants_bck(self, x1, x2, sequence, polymorphismsKeys) :
		ret = []
		for pk in polymorphismsKeys :
			if pk >= x1 or pk < x2 :
				for c in self.polymorphisms[pk] :
					sequence[pk-x1] = c
					ret.extend(self.__getSequenceVariants(self, x1, x2, sequence, polymorphismsKeys))
	
	def __getSequenceVariants(self, x1, polyStart, polyStop, listSequence) :
		"""polyStop, is the polymorphisme number where te calcul of combinaisons stops"""
		if polyStart < len(self.polymorphisms) and polyStart < polyStop: 
			sequence = copy.copy(listSequence)
			ret = []
			
			polyStartNext = polyStart
			pk = self.polymorphisms[polyStartNext]
			if pk[0]-x1 < len(listSequence) : 
				while pk[0] < x1 :
					pk = self.polymorphisms[polyStartNext]
					polyStartNext += 1
					
				ret.extend(self.__getSequenceVariants(x1, polyStartNext +1, polyStop, sequence))
				for allele in pk[1][:-1] :
					sequence[pk[0]-x1] = allele
					ret.extend(self.__getSequenceVariants(x1, polyStartNext +1, polyStop, sequence))
			else :
				return [''.join(listSequence)]
		else :
			return [''.join(listSequence)]
			
		return ret
	
	def getSequenceVariants(self, x1 = 0, x2 = -1, maxVariantNumber = 128) :
		"""compute the combinaisons of all polymorphismes between x1 and x2
		@return a couple (bool, variants of sequence between x1 and x2), bool == true if the maxVariantNumber threshold is depasse"""
		if x2 == -1 :
			xx2 = len(self.defaultSequence)
		else :
			xx2 = x2
		
		polyStop = 0
		nbP = 1
		stopped = False
		for p in self.polymorphisms:
			if x1 <= p[0] and p[0] <= xx2 :
				nbP *= len(p[1])
				if nbP > maxVariantNumber :
					stopped = True
					break
				polyStop += 1
		#print polyStop, nbP
		#if nbP <= maxVariantNumber :
		return (stopped, self.__getSequenceVariants(x1, 0, polyStop, list(self.defaultSequence[x1:xx2])))
	
	
	def getNbVariants(self, x1, x2 = -1) :
		"""@returns the nb of variants of sequences between x1 and x2"""
		if x2 == -1 :
			xx2 = len(self.defaultSequence)
		else :
			xx2 = x2
		
		nbP = 1
		for p in self.polymorphisms:
			if x1 <= p[0] and p[0] <= xx2 :
				nbP *= len(p[1])
		
		return nbP
		
	def find(self, strSeq) :
		lpos = range(len(self))
		arr = self.encode(strSeq)
		
		if len(arr) == 0 :
			print "strSeq encodeion is empty", strSeq, arr
			return []

		for i in range(len(arr)) :
			l = []
			for j in range(len(lpos)) :
				#print i, arr[i], self[lpos[j]+i], arr[i] & self[lpos[j]+i]
				if (lpos[j]+i < len(self)) and (arr[i] & self[lpos[j]+i] > 0):
					l.append(lpos[j])
			lpos = l
		
		return l
		
	def find_naive(self, strSeq) :
		
		res = -1
		arr = self.encode(strSeq)
		
		#naive search, faire qlq de plus fancy quand j'aurais le temps	
		for i in range(len(self)) :
			found = True
			if (i+len(arr) <= len(self)) :
				for j in range(len(arr)) :
					if (self[i+j] & arr[j]) == 0 :
						found = False
						break
				if found :
					return i
				
		return -1
		
	def findAll_naive(self, strSeq) :
	
		arr = self.encode(strSeq)
		self.__testFind(arr)
		#naive search, faire qlq de plus fancy quand j'aurais le temps
		ret = []
		for i in range(len(self)) :
			found = True
			if (i+len(arr) <= len(self)) :
				for j in range(len(arr)) :
					if (self[i+j] & arr[j]) == 0 :
						found = False
						break
				if found :
					ret.append(i)
		
		return ret
	
	def __and__(self, arr) :
		self.__testBinary(arr)
		
		ret = BinarySequence(self.typecode, self.forma, self.charToBin)
		for i in range(len(arr)) :
			ret.append(self[i] & arr[i])
		
		return ret
	
	def __xor__(self, arr) :
		self.__testBinary(arr)
		
		ret = BinarySequence(self.typecode, self.forma, self.charToBin)
		for i in range(len(arr)) :
			ret.append(self[i] ^ arr[i])
		
		return ret

	def __eq__(self, seq) :
		self.__testBinary(seq)
		
		if len(seq) != len(self) :
			return False

		return all( self[i] == seq[i] for i in range(len(self)) )

		
	def append(self, arr) :
		self.binSequence.append(arr)

	def extend(self, arr) :
		self.binSequence.extend(arr)

	def decode(self, binSeq):
		ret = ''
		for b in binSeq :
			ch = ''
			for c in self.charToBin :
				if b & self.forma[self.charToBin[c]] > 0 :
					ch += c +'/'
			if ch == '' :
				raise KeyError('Key %d unkowom, bad format' % b)
			
			ret += ch[:-1] #remove last '/'
		return ret
		
	def getChar(self, i):
		return self.decode([self.binSequence[i]])
		
	def __len__(self):
		return len(self.binSequence)

	def __getitem__(self,i):
		return self.binSequence[i]
	
	def __setitem__(self, i, v):
		self.binSequence[i] = v

class AABinarySequence(BinarySequence) :
	def __init__(self, sequence):
		f = array.array('I', [1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L, 256L, 512L, 1024L, 2048L, 4096L, 8192L, 16384L, 32768L, 65536L, 131072L, 262144L, 524288L, 1048576L])
		c = {'A': 17, 'C': 14, 'E': 19, 'D': 15, 'G': 13, 'F': 16, 'I': 3, 'H': 9, 'K': 8, '*': 1, 'M': 20, 'L': 0, 'N': 4, 'Q': 11, 'P': 6, 'S': 7, 'R': 5, 'T': 2, 'W': 10, 'V': 18, 'Y': 12}
		BinarySequence.__init__(self, sequence, f, c)
	
class NucBinarySequence(BinarySequence) :
	def __init__(self, sequence):
		f = array.array('B', [1, 2, 4, 8])
		c = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
		self.nucleotides = {
			'A' : ['A'], 'T' : ['T'], 'C' : ['C'], 'G' : ['G'],
			'R' : ['A','G'], 'Y' : ['C','T'], 'M': ['A','C'],
			'K' : ['T','G'], 'W' : ['A','T'], 'S' : ['C','G'], 'B': ['C','G','T'],
			'D' : ['A','G','T'], 'H' : ['A','C','T'], 'V' : ['A','C','G'], 'N': ['A','C','G','T']
			}
		BinarySequence.__init__(self, sequence, f, c)
		
	def encode(self, sequence):
		res = array.array(self.forma.typecode)
		i = 0
		for c in sequence :
			b = 0
			for nuc in self.nucleotides[c] :
				b = b | self.forma[self.charToBin[nuc]]
			if c in uf.polymorphicNucleotides :
				self.polymorphisms.append((i, self.nucleotides[c]))
			
			res.append(b)
			i += 1
		return res
		
if __name__=="__main__":
	a = AABinarySequence('ACEDGFIHKM')
	print a.find_fast('ED')
