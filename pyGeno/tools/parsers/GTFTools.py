import gzip

class GTFEntry(object) :
	def __init__(self, gtfFile, lineNumber) :
		"""A single entry in a GTF file"""
		
		self.lineNumber = lineNumber
		self.gtfFile = gtfFile
		self.data = gtfFile.lines[lineNumber][:-2].split('\t') #-2 remove ';\n'
		proto_atts = self.data[gtfFile.legend['attributes']].strip().split('; ')
		atts = {}
		for a in proto_atts :
			sa = a.split(' ')
			atts[sa[0]] = sa[1].replace('"', '')
		self.data[gtfFile.legend['attributes']] = atts
	
	def __getitem__(self, k) :
		try :
			return self.data[self.gtfFile.legend[k]]
		except KeyError :
			try :
				return self.data[self.gtfFile.legend['attributes']][k]
			except KeyError :
				#return None
				raise KeyError("Line %d does not have an element %s.\nline:%s" %(self.lineNumber, k, self.gtfFile.lines[self.lineNumber]))
	
	def __repr__(self) :
		return "<GTFEntry line: %d>" % self.lineNumber
	
	def __str__(self) :
		return  "<GTFEntry line: %d, %s>" % (self.lineNumber, str(self.data))

class GTFFile(object) :
	"""This is a simple GTF2.2 (Revised Ensembl GTF) parser, see http://mblab.wustl.edu/GTF22.html for more infos"""
	def __init__(self, filename, gziped = False) :
		
		self.filename = filename
		self.legend = {'seqname' : 0, 'source' : 1, 'feature' : 2, 'start' : 3, 'end' : 4, 'score' : 5, 'strand' : 6, 'frame' : 7, 'attributes' : 8}

		if gziped : 
			f = gzip.open(filename)
		else :
			f = open(filename)
		
		self.lines = []
		for l in f :
			if l[0] != '#' and l != '' :
				self.lines.append(l)
		f.close()
		
		self.currentIt = -1

	def get(self, line, elmt) :
		"""returns the value of the field 'elmt' of line 'line'"""
		return self[line][elmt]

	def __iter__(self) :
		self.currentPos = -1
		return self

	def next(self) :
		self.currentPos += 1
		try :
			return GTFEntry(self, self.currentPos)
		except IndexError:
			raise StopIteration

	def __getitem__(self, i) :
		"""returns the ith entry"""
		if self.lines[i].__class__ is not GTFEntry :
			self.lines[i] = GTFEntry(self, i)
		return self.lines[i]

	def __repr__(self) :
		return "<GTFFile: %s>" % (os.path.basename(self.filename))

	def __str__(self) :
		return "<GTFFile: %s, gziped: %s, len: %d>" % (os.path.basename(self.filename), self.gziped, len(self))
	
	def __len__(self) :
		return len(self.lines)
