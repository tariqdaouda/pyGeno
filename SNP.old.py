import tools

class SNPVariant :
	#'gene-symbol;fxn-class;allele;residue;aa_position;frame'
	def __init__(self, snp, geneSymbol, function, allele, residue, aaPosition, frame) :
		self.geneSymbol = geneSymbol
		self.snp = snp
		self.function = function
		self.allele = allele
		self.residue = residue
		
		if aaPosition == '' :
			self.aaPosition = 0
		else :
			self.aaPosition = int(aaPosition)
		
		self.frame = frame
		
		self.eqTest = {}
		self.eqTest[self.toCSV()] = 0
		
	def __eq__(self, otherVariant):
		"""returns True if otherVariant has the same values, False if not"""
		try :
			self.eqTest[otherVariant.toCSV()] += 1
			return True
		except :
			return False
		
	def toXML(self) :
		return '<variant geneSymbol="%s" function="%s" allele="%s" residue="%s" aaPosition="%s" frame="%s"/>'%(self.geneSymbol, self.function, self.allele, self.residue, self.aaPosition, self.frame)
	
	def toCSV(self) :
		return '%s;%s;%s;%s;%s;%s' % (self.geneSymbol, self.function, self.allele, self.residue, self.aaPosition, self.frame)

class SNP :
	"""A SNP"""
	#'chr;chr-pos;rsid;alleles;orient;assembly;validated'
	def __init__(self, chromosome, x1, rsId, alleles, orient, assembly, validated) :
		"""A snp from dbSNP"""
		
		self.chromosome = chromosome
		self.x1 = int(x1)
		self.rsId = rsId
		self.alleles = alleles
		self.orient = orient
		self.assembly = assembly
		self.validated = (validated.strip().lower() == 'yes')
		self.variants = []
		
	def addVariant(self, geneSymbol, function, allele, residue, aaPosition, frame) :
		"""Adds variant to the SNP if it hasn't already been added"""
		v = SNPVariant(self, geneSymbol, function, allele, residue, aaPosition, frame)
		if not self.hasVariant(v) :
			self.variants.append(v)
			#print "add"
			
	def hasVariant(self, variant) :
		"""Returns True if the varians is aready stored, False if not"""
		for v in self.variants :
			if v == variant :
				return True
		return False
		
	def __getitem__(self, i) :
		"""returns the ith variant"""
		return self.variants[i]
	
	def getVariants(self, function = '') :
		"""function is the ncbi function {missense, stop-gained, 
		reference, intron-variant, synonymous-codon,  frameshift-variant,
		utr-variant-3-prime, utr-variant-5-prime, downstream-variant-*, upstream-variant-*,
		etc...}"""
		if function == '' :
			return self.variants
		else :
			vars = []
			for v in self.variants :
				if v.function == function :
					vars.append[v]
			return vars
	
	def getReferenceVariant(self) :
		ret = None
		for v in self.variants :
			if v.function == 'reference' :
				return v
	
	def getCsvLegend(self, separator) :
		pass

	def CSVForVariant(self, i):
		return '"%s";"%s";"%s";"%s";"%s";"%s";"%s";%s'%(self.chromosome, self.x1, self.rsId, self.alleles, self.orient, self.assembly, self.validated, self.variants[i].toCSV())
		
	def toCSV(self) :
		"""One line per variant"""
		csv = ''
		for i in range(len(self.variants)) :
			csv += self.CSVForVariant(i) + '\n'
		return csv
		
	def _getVariantsXMl(self) :
		str = ''
		for v in self.variants :																
			str += v.toXML()
			
		return str
	
	def toXML(self) :
		""""returns a proper xml representation of the object """
		variants = ""
		for v in self.variants :
			variants += v.toXML()
		
		#'chr;chr-pos;rsid;alleles;orient;assembly;validated'
		return '<snp chromosome="%s" x1="%s" rsId="%s" alleles="%s" orient="%s" assembly="%s" validated="%s">%s</snp>'%(self.chromosome, self.x1, self.rsId, self.alleles, self.orient, self.assembly, self.validated, variants)
	
	def __str__(self) :
		str = self.toXML()
		str = str.replace('<snp>', '').replace('</snp>', '').replace('<variant', '\n\t> variant').replace('<', '').replace('/>', '')
		return str

