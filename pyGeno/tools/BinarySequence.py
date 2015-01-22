import array, copy
import UsefulFunctions as uf

class BinarySequence :
	"""A class for representing sequences in a binary format"""
	
	def __init__(self, sequence, arrayForma, charToBinDict) :
	
		self.forma = arrayForma
		self.charToBin = charToBinDict
		self.sequence = sequence
		
		self.binSequence, self.defaultSequence, self.polymorphisms = self.encode(sequence)
		self.itemsize = self.binSequence.itemsize
		self.typecode = self.binSequence.typecode
		#print 'bin', len(self.sequence), len(self.binSequence)

	def encode(self, sequence):
		"""Returns a tuple (binary reprensentation, default sequence, polymorphisms list)"""
		
		polymorphisms = []
		defaultSequence = '' 
		binSequence = array.array(self.forma.typecode)
		b = 0
		i = 0
		trueI = 0 #not inc in case if poly
		poly = set()
		while i < len(sequence)-1:
			b = b | self.forma[self.charToBin[sequence[i]]]
			if sequence[i+1] == '/' :
				poly.add(sequence[i])
				i += 2
			else :
				binSequence.append(b)
				if len(poly) > 0 :
					poly.add(sequence[i])
					polymorphisms.append((trueI, poly))
					poly = set()
				
				bb = 0
				while b % 2 != 0 :
					b = b/2
					
				defaultSequence += sequence[i]
				b = 0
				i += 1
				trueI += 1
		
		if i < len(sequence) :
			b = b | self.forma[self.charToBin[sequence[i]]]
			binSequence.append(b)
			if len(poly) > 0 :
				if sequence[i] not in poly :
					poly.add(sequence[i])
				polymorphisms.append((trueI, poly))
			defaultSequence += sequence[i]
		
		return (binSequence, defaultSequence, polymorphisms)

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
		"""
		Compares strSeq with self.sequence.
		If not 'strict', this function ignores the cases of matching heterozygocity (ex: for a given position i, strSeq[i] = A and self.sequence[i] = 'A/G'). If 'strict' it returns all positions where strSeq differs self,sequence
		"""
		arr = self.encode(strSeq)[0]
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
		"""returns all polymorphsims in the form of a dict pos => alleles"""
		return self.polymorphisms
	
	def getDefaultSequence(self) :
		"""returns a default version of sequence where only the last allele of each polymorphisms is shown"""
		return self.defaultSequence
	
	def __getSequenceVariants(self, x1, polyStart, polyStop, listSequence) :
		"""polyStop, is the polymorphisme at wixh number where the calcul of combinaisons stops"""
		if polyStart < len(self.polymorphisms) and polyStart < polyStop: 
			sequence = copy.copy(listSequence)
			ret = []
			
			pk = self.polymorphisms[polyStart]
			posInSequence = pk[0]-x1
			if posInSequence < len(listSequence) : 
				for allele in pk[1] :
					sequence[posInSequence] = allele
					ret.extend(self.__getSequenceVariants(x1, polyStart+1, polyStop, sequence))
			
			return ret
		else :
			return [''.join(listSequence)]

	def getSequenceVariants(self, x1 = 0, x2 = -1, maxVariantNumber = 128) :
		"""returns the sequences resulting from all combinaisons of all polymorphismes between x1 and x2. The results is a couple (bool, variants of sequence between x1 and x2), bool == true if there's more combinaisons than maxVariantNumber"""
		if x2 == -1 :
			xx2 = len(self.defaultSequence)
		else :
			xx2 = x2
		
		polyStart = None
		nbP = 1
		stopped = False
		j = 0
		for p in self.polymorphisms :
			if p[0] >= xx2 :
				break
			
			if x1 <= p[0] :
				if polyStart == None :
					polyStart = j
				
				nbP *= len(p[1])
				
				if nbP > maxVariantNumber :
					stopped = True
					break
			j += 1
		
		if polyStart == None :
			return (stopped, [self.defaultSequence[x1:xx2]])
		
		return (stopped, self.__getSequenceVariants(x1, polyStart, j, list(self.defaultSequence[x1:xx2])))

	def getNbVariants(self, x1, x2 = -1) :
		"""returns the nb of variants of sequences between x1 and x2"""
		if x2 == -1 :
			xx2 = len(self.defaultSequence)
		else :
			xx2 = x2
		
		nbP = 1
		for p in self.polymorphisms:
			if x1 <= p[0] and p[0] <= xx2 :
				nbP *= len(p[1])
		
		return nbP

	def _dichFind(self, needle, currHaystack, offset, lst = None) :
		"""dichotomic search, if lst is None, will return the first position found. If it's a list, will return a list of all positions in lst. returns -1 or [] if no match found"""
		
		if len(currHaystack) == 1 :
			if (offset <= (len(self) - len(needle))) and (currHaystack[0] & needle[0]) > 0 and (self[offset+len(needle)-1] & needle[-1]) > 0 :
				found = True
				for i in xrange(1, len(needle)-1) :
					if self[offset + i] & needle[i] == 0 :
						found = False
						break
				if found :
					if lst is not None :
						lst.append(offset)
					else :
						return offset
				else :
					if lst is None :
						return -1
		else :
			if (offset <= (len(self) - len(needle))) :
				if lst is not None :
					self._dichFind(needle, currHaystack[:len(currHaystack)/2], offset, lst)
					self._dichFind(needle, currHaystack[len(currHaystack)/2:], offset + len(currHaystack)/2, lst)
				else :
					v1 = self._dichFind(needle, currHaystack[:len(currHaystack)/2], offset, lst)
					if v1 > -1 :
						return v1
					
					return self._dichFind(needle, currHaystack[len(currHaystack)/2:], offset + len(currHaystack)/2, lst)
			return -1

	def find(self, strSeq) :
		"""returns the first occurence of strSeq in self. Takes polymorphisms into account"""
		arr = self.encode(strSeq)
		return self._dichFind(arr[0], self, 0, lst = None)

	def findAll(self, strSeq) :
		"""Same as find but returns a list of all occurences"""
		arr = self.encode(strSeq)
		lst = []
		self._dichFind(arr[0], self, 0, lst)
		return lst
		
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

	def decode(self, binSequence):
		"""decodes a binary sequence to return a string"""
		try:
			binSeq = iter(binSequence[0])
		except TypeError, te:
			binSeq = binSequence
    
		ret = ''
		for b in binSeq :
			ch = ''
			for c in self.charToBin :
				if b & self.forma[self.charToBin[c]] > 0 :
					ch += c +'/'
			if ch == '' :
				raise KeyError('Key %d unkowom, bad format' % b)
			
			ret += ch[:-1]
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
	"""A binary sequence of amino acids"""
	
	def __init__(self, sequence):
		f = array.array('I', [1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L, 256L, 512L, 1024L, 2048L, 4096L, 8192L, 16384L, 32768L, 65536L, 131072L, 262144L, 524288L, 1048576L])
		c = {'A': 17, 'C': 14, 'E': 19, 'D': 15, 'G': 13, 'F': 16, 'I': 3, 'H': 9, 'K': 8, '*': 1, 'M': 20, 'L': 0, 'N': 4, 'Q': 11, 'P': 6, 'S': 7, 'R': 5, 'T': 2, 'W': 10, 'V': 18, 'Y': 12}
		BinarySequence.__init__(self, sequence, f, c)
	
class NucBinarySequence(BinarySequence) :
	"""A binary sequence of nucleotides"""
	
	def __init__(self, sequence):
		f = array.array('B', [1, 2, 4, 8])
		c = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
		ce = {
			'R' : 'A/G', 'Y' : 'C/T', 'M': 'A/C',
			'K' : 'T/G', 'W' : 'A/T', 'S' : 'C/G', 'B': 'C/G/T',
			'D' : 'A/G/T', 'H' : 'A/C/T', 'V' : 'A/C/G', 'N': 'A/C/G/T'
			}
		lstSeq = list(sequence)
		for i in xrange(len(lstSeq)) :
			if lstSeq[i] in ce :
				lstSeq[i] = ce[lstSeq[i]]
		lstSeq = ''.join(lstSeq)
		BinarySequence.__init__(self, lstSeq, f, c)

if __name__=="__main__":
	
	def test0() :
		#seq = 'ACED/E/GFIHK/MLMQPS/RTWVY'
		seq = 'ACED/E/GFIHK/MLMQPS/RTWVY/A/R'
		bSeq = AABinarySequence(seq)
		start = 0
		stop = 4
		rB = bSeq.getSequenceVariants_bck(start, stop)
		r = bSeq.getSequenceVariants(start, stop)
		
		#print start, stop, 'nb_comb_r', len(r[1]), set(rB[1])==set(r[1]) 
		print start, stop#, 'nb_comb_r', len(r[1]), set(rB[1])==set(r[1]) 
		
		#if set(rB[1])!=set(r[1]) :
		print '-AV-'
		print start, stop, 'nb_comb_r', len(rB[1])
		print '\n'.join(rB[1])
		print '=AP========'
		print start, stop, 'nb_comb_r', len(r[1]) 
		print '\n'.join(r[1])
	
	def testVariants() :
		seq = 'ATGAGTTTGCCGCGCN'
		bSeq = NucBinarySequence(seq)
		print bSeq.getSequenceVariants() 

	testVariants()
