import tools
import copy

class InvalidProtein(Exception):
	def __init__(self, peptide, proteinId, protein, gene, message):
		self.peptide = peptide
		self.proteinId = proteinId
		self.protein = protein
		self.gene = gene
		self.message = message
	
	def __str__(self):
		return """
		Description : %s
		peptide : %s
		protein_id : %s
		protein : %s
		gene_symbol : %s"""%(self.message, self.peptide, self.proteinId, self.protein, self.gene.symbol)
		
class RefactorizedExons :
	""""A data structure that represent a refactorized exon list, in other terms a list of exons where the
	the length of each exons is a multiple of 3 """
	def __init__(self, gene, protId):
		self.gene = gene
		self.protein = protId
		
		self.exonSequences = self.gene.getCdsSequences(protId)
		self.exonPosition = self.gene.getCdsPosistions(protId)
		self._refactorizeExons()
		
	def _refactorizeExons(self) :
	"""Create refact versions of the exon sequence and positions, by rearging them according to the gene direction
	and making the lenght of all exon multiples of 3"""
	
		self.refactExonSequences = copy.copy(self.exonSequences)
		self.refactExonPositions = copy.copy(self.exonPositions)
		
		self.refactExonTranslatedSequences = copy.copy(self.exonSequences)
		
		#Direction
		if self.gene.direction == '-' :
			for i in range(len(self.refactExonSequences)) :
				self.refactExonSequences[i] = tools.complement(self.exonSequences[i])
				self.refactExonPositions[i] = self.exonPositions[i][::-1]
			
			protExon = self.refactExonSequences[::-1]
			protPos = self.refactExonPositions[::-1]
			
		#multiples of 3		
		for i in range(len(self.refactExonSequences)) :
			mod = len(self.refactExonSequences[i])%3
			
			if mod != 0 :
				exc = self.refactExonSequences[i][-mod:]
				self.refactExonSequences[i] = self.refactExonSequences[i][:-mod]
	
				if self.direction == '-' :
					self.refactExonPositions[i] = (self.refactExonPositions[i][0], self.refactExonPositions[i][1] + mod)
				elif self.direction == '+' :
					self.refactExonPositions[i] = (self.refactExonPositions[i][0], self.refactExonPositions[i][1] - mod)
				
				if i < len(self.refactExonSequences)-1 :
					self.refactExonSequences[i+1]  = exc + self.refactExonSequences[i+1]
					
			self.refactExonTranslatedSequences[i] = tools.translateDNA(self.refactExonSequences[i])

		def __getitem__(self, i)
			return self.refactorizedExon[i]
		
class ExtendedList :
	"""This data structure augmentes a list by creating element between elements that are combination of the
	suffix of the previous element and prefix of the following one"""
	def __init__(self, splitSequence, overlapRange) :
		"""splitSequence : the list to augment
		overlapRange : the length of the suffix/prefix"""
		
		self.overlapRange = overlapRange
		
		lenSeq = len(splitSequence)
		self.tree = range(lenSeq*2-1)
		for i in range(lenSeq) :
			self.tree[i*2] = splitSequence[i]
		
		for i in range(1, len(self.tree) -1, 2) :
			self.tree[i] = self.tree[i-1][-overlapRange:] + self.tree[i+1][:overlapRange]
	
	def __getitem__(self, i) :
		return self.tree[i]
	
	def __setitem__(self, i, v) :
		self.tree[i] = v

	def __len__(self) :
		return len(self.tree)

	def __str__(self) :
		str = ''
		for i in range(len(self.tree)) :
			if i%2 == 0 :
				str += '==> |'+self.tree[i]+'|\n'
			if i%2 ==  1 :
				str += '-----> |'+self.tree[i]+'|\n'
		
		return str
		
class PositionWizard :

	def __init__(self) :
		
		self.gene = gene
		self.proteinId = proteinId
		
		self.cdsSequences = gene.getProteinCdsSequences(proteinId)
		self.cdsPositions = gene.getProteinCdsPosistions(proteinId)
		self.direction = gene.direction
		
		self._refactorizeCds()
		self.protein = ''.join(self.refactCdsTranslatedSequences)
		

	def proteinToCDNA(self, pos) :
		return pos*3
	
	def cDNAtoProtein(self, pos) :
		return int(pos/3)
		
	def proteinToDNA(self, pos, protein) :
		refact = RefactorizedExons(protein.gene, protein.id)
		
		tree = ExtendedList(self.refactCdsTranslatedSequences, len(peptide) - 1)
		
		positions = []
		
		for i in range(len(tree)) :
			
			nucPepRatio = 3
			pepPositions = tools.findAll(tree[i], peptide)
			for pepPos in pepPositions :

				if pepPos >= 0 :
					if i%2 == 0 :				
						if (self.direction == '+') :
							nucPlus = len(self.refactCdsSequences[i/2]) - len(self.cdsSequences[i/2])
							
							dnaStartPos = pepPos*nucPepRatio + self.refactCdsPositions[i/2][0]
							if  nucPlus > 0 and pepPos == 0: #si le pep est au debut, une partie fait parti du cds precedant
								positions.append((self.cdsPositions[i/2 - 1][1] - nucPlus - 1, self.cdsPositions[i/2 - 1][1]))
								dnaStopPos = dnaStartPos+len(peptide)*nucPepRatio-1 - nucPlus #On enleve le nombre de nuc apparetant au cds precedant
							else :
								dnaStopPos = dnaStartPos+len(peptide)*nucPepRatio-1
							
						elif (self.direction == '-') :
							nucPlus = len(self.refactCdsSequences[i/2]) - len(self.cdsSequences[i/2])
							
							dnaStartPos = self.refactCdsPositions[i/2][0] - pepPos*nucPepRatio
							if  nucPlus > 0 and (pepPos + len(peptide) - 1) == len(tree[i]):
								positions.append((self.cdsPositions[i/2 + 1][0], self.cdsPositions[i/2 + 1][0] + nucPlus - 1))
								dnaStopPos = dnaStartPos-len(peptide)*nucPepRatio+1 + nucPlus
							else :
								dnaStopPos = dnaStartPos-len(peptide)*nucPepRatio+1	
							
							
						positions.append((dnaStartPos, dnaStopPos))
						
					else :
						
						nucPlus = len(self.refactCdsSequences[i/2]) - len(self.cdsSequences[i/2])
	
						if (self.direction == '+') :
							dnaStartPos = self.refactCdsPositions[i/2][1] - ((len(tree[i])/2)-pepPos)*nucPepRatio + 1 + nucPlus
							dnaStopPos = self.refactCdsPositions[i/2+1][0] + (len(peptide) - (len(tree[i])/2 - pepPos))*nucPepRatio -1 - nucPlus
							
							positions.append((dnaStartPos, self.refactCdsPositions[i/2][1]))
							positions.append((self.refactCdsPositions[i/2+1][0], dnaStopPos))
							
						else : 
							nucStartPos = self.refactCdsPositions[i/2][1] + ((len(tree[i])/2)-pepPos)*nucPepRatio - 1 + nucPlus
							nucStopPos = self.refactCdsPositions[i/2+1][0] - (len(peptide) - (len(tree[i])/2 - pepPos))*nucPepRatio - nucPlus
													
							positions.append((self.refactCdsPositions[i/2][1], nucStartPos))
							positions.append((nucStopPos, self.refactCdsPositions[i/2+1][0]))
		
		return positions

	def test(self) :
		prot = ''.join(self.refactCdsTranslatedSequences)
		pepLen = 10

		print "====\nTesting with protein: %s\n\n"%(prot)
		for pepPos in range(0, len(prot), 1) :
			print "----->positions:", self.refactCdsPositions
			peptide = prot[pepPos:pepPos + pepLen]
			print "--->peptide '%s'"%(peptide), self.peptidFindToNucleotidPos(peptide, True), '\n=====\n'
