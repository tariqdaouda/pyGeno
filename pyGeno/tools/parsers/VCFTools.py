import os, types, gzip

class VCFEntry(object) :
	"""A single entry in a VCF file"""
	
	def __init__(self, vcfFile, line, lineNumber) :
		#CHROM POS     ID        REF    ALT     QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		#20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
		self.vcfFile = vcfFile
		self.lineNumber = lineNumber
		self.data = {}
		
		tmpL = line.replace('\r', '\n').replace('\n', '')
		tmpData = str(tmpL).split('\t')
		for i in range(6) :
			self.data[vcfFile.dnegel[i]] = tmpData[i]
		self.data['POS'] = int(self.data['POS'])
		
		filters = tmpData[6].split(';')
		if len(filters) == 1 :
			self.data['FILTER'] = filters
		else :
			for filter_value in tmpData[6].split(';') :
				filt, value = info_value.split('=')
				self.data['FILTER'][filt] = value
		
		self.data['INFO'] = {}
		for s in tmpData[7].split(';') :
			info_value = s.split('=')
			try :
				typ = self.vcfFile.meta['INFO'][info_value[0]]['Type']
			except KeyError :
				typ = None
				
			if len(info_value) == 1 :
				if typ == 'Flag' or typ == None :
					self.data['INFO'][info_value[0]] = True
				else :
					raise ValueError('%s is not a flag and has no value' % info_value[0])
			else :
				if typ == 'Integer' :
					self.data['INFO'][info_value[0]] = int(info_value[1])
				elif typ == 'Float' :
					self.data['INFO'][info_value[0]] = float(info_value[1])
				else :
					self.data['INFO'][info_value[0]] = info_value[1]
				
	def __getitem__(self, key) :
		"with the vcf file format some fields are not present in all elements therefor, this fct never raises an exception but returns None or False if the field is definied as a Flag in Meta"
		try :
			return self.data[key]
		except KeyError:
			try :
				return self.data['INFO'][key]
			except KeyError:
				try :
					if self.vcfFile.meta['INFO'][key]['Type'] == 'Flag' :
						self.data['INFO'][key] = False
						return self.data['INFO'][key]
					else :
						return None
				except KeyError:
					return None
	
	def __repr__(self) :
		return "<VCFEntry line: %d>" % self.lineNumber
	
	def __str__(self) :
		return "<VCFEntry line: %d,  %s" % (self.lineNumber, str(self.data))

class VCFFile(object) :
	"""
	This is a small parser for VCF files, it should work with any VCF file but has only been tested on dbSNP138 files.
	Represents a whole VCF file::
		
		#reading
		f = VCFFile()
		f.parse('hop.vcf')
		for line in f :
			print line['pos']
	"""
	
	def __init__(self, filename = None, gziped = False, stream = False) :
		self.legend = {}
		self.dnegel = {}
		self.meta = {}
		self.lines = None
		if filename :
			self.parse(filename, gziped, stream)
		
	def parse(self, filename, gziped = False, stream = False) :
		"""opens a file"""
		self.stream = stream
		
		if gziped :
			self.f = gzip.open(filename)
		else :
			self.f = open(filename)
		
		self.filename = filename
		self.gziped = gziped
		
		lineId = 0
		inLegend = True
		while inLegend :
			ll = self.f.readline()
			l = ll.replace('\r', '\n').replace('\n', '')
			if l[:2] == '##' :
				eqPos = l.find('=')
				key = l[2:eqPos]
				values = l[eqPos+1:].strip()
				
				if l[eqPos+1] != '<' :
					self.meta[key] = values
				else :
					if key not in self.meta :
						self.meta[key] = {}
					svalues = l[eqPos+2:-1].split(',') #remove the < and > that surounf the string 
					idKey = svalues[0].split('=')[1]
					self.meta[key][idKey] = {} 
					i = 1
					for v in svalues[1:] :
						sv = v.split("=")
						field = sv[0]
						value = sv[1]
						if field.lower() == 'description' :
							self.meta[key][idKey][field] = ','.join(svalues[i:])[len(field)+2:-1]
							break
						self.meta[key][idKey][field] = value
						i += 1
			elif l[:6] == '#CHROM': #we are in legend
				sl = l.split('\t')
				for i in range(len(sl)) :
					self.legend[sl[i]] = i
					self.dnegel[i] = sl[i]
				break
			
			lineId += 1
		
		if not stream :
			self.lines = self.f.readlines()
			self.f.close()
	
	def close(self) :
		"""closes the file"""
		self.f.close()
		
	def _developLine(self, lineNumber) :
		if self.lines[lineNumber].__class__ is not VCFEntry :
			self.lines[lineNumber] = VCFEntry(self, self.lines[lineNumber], lineNumber)
		
	def __iter__(self) :
		self.currentPos = -1
		return self
	
	def next(self) :
		self.currentPos += 1
		if not self.stream :
			try :
				return self[self.currentPos-1]
			except IndexError:
				raise StopIteration
		else :
			line = self.f.readline()
			if not line :
				raise StopIteration
			return VCFEntry(self, line, self.currentPos)
	
	def __getitem__(self, line) :
		"""returns the lineth element"""
		if self.stream :
			raise KeyError("When the file is opened as a stream it's impossible to ask for specific item")
		
		if self.lines[line].__class__ is not VCFEntry :
			self._developLine(line)
		return self.lines[line]

	def __len__(self) :
		"""returns the number of entries"""
		return len(self.lines)
	
	def __repr__(self) :
		return "<VCFFile: %s>" % (os.path.basename(self.filename))
	
	def __str__(self) :
		if self.stream :
			return "<VCFFile: %s, gziped: %s, stream: %s, len: undef>" % (os.path.basename(self.filename), self.gziped, self.stream)
		else :
			return "<VCFFile: %s, gziped: %s, stream: %s, len: %d>" % (os.path.basename(self.filename), self.gziped, self.stream, len(self))
		
if __name__ == '__main__' :
	from pyGeno.tools.ProgressBar import ProgressBar
	
	#v = VCFFile('/u/daoudat/py/pySrc/pyGeno_stuff/V2/packages/dbSNP/human/dbSNP138/00-All.vcf.gz', gziped = True, stream = True)
	v = VCFFile('/u/daoudat/py/pySrc/pyGeno_stuff/V2/packages/dbSNP/human/dbSNP138/test/00-All-test.vcf.gz', gziped = True, stream = True)
	startTime = time.time()
	i = 0
	pBar = ProgressBar()
	for f in v :
		#print f
		pBar.update('%s' % i)
		if i > 1000000 :
			break
		i += 1
	pBar.close()
	
	#v = VCFFile('/u/daoudat/py/pySrc/pyGeno_stuff/V2/packages/dbSNP/human/dbSNP138/00-All.vcf', gziped = False, stream = True)
	v = VCFFile('/u/daoudat/py/pySrc/pyGeno_stuff/V2/packages/dbSNP/human/dbSNP138/test/00-All-test.vcf', gziped = False, stream = False)
	startTime = time.time()
	i = 0
	pBar = ProgressBar()
	for f in v :
		pBar.update('%s' % i)
		if i > 1000000 :
			break
		i += 1
		#print f
	pBar.close()
	#print v.lines

