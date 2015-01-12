import string, os, copy, types

class UnknownNucleotide(Exception) :
	def __init__(self, nuc) :
		self.msg =  'Unknown nucleotides %s' % str(nuc)

	def __str__(self) :
		return self.msg

#This will probably be moved somewhere else in the futur
def saveResults(directoryName, fileName, strResults, log = '', args = ''):

	if not os.path.exists(directoryName):
		os.makedirs(directoryName)

	resPath = "%s/%s"%(directoryName, fileName)
	resFile = open(resPath, 'w')
	print "Saving results :\n\t%s..."%resPath
	resFile.write(strResults)
	resFile.close()

	if log != '' :
		errPath = "%s.err.txt"%(resPath)
		errFile = open(errPath, 'w')

		print "Saving log :\n\t%s..." %errPath
		errFile.write(log)
		errFile.close()

	if args != '' :
		paramPath = "%s.args.txt"%(resPath)
		paramFile = open(paramPath, 'w')

		print "Saving arguments :\n\t%s..." %paramPath
		paramFile.write(args)
		paramFile.close()

	return "%s/"%(directoryName)

nucleotides = ['A', 'T', 'C', 'G']
polymorphicNucleotides = {
			'R' : ['A','G'], 'Y' : ['C','T'], 'M': ['A','C'],
			'K' : ['T','G'], 'W' : ['A','T'], 'S' : ['C','G'], 'B': ['C','G','T'],
			'D' : ['A','G','T'], 'H' : ['A','C','T'], 'V' : ['A','C','G'], 'N': ['A','C','G','T']
			}

#<7iyed>
#from Molecular Systems Biology 8; Article number 572; doi:10.1038/msb.2012.3
codonAffinity = {'CTT': 'low', 'ACC': 'high', 'ACA': 'low', 'ACG': 'high', 'ATC': 'high', 'AAC': 'high', 'ATA': 'low', 'AGG': 'high', 'CCT': 'low', 'ACT': 'low', 'AGC': 'high', 'AAG': 'high', 'AGA': 'low', 'CAT': 'low', 'AAT': 'low', 'ATT': 'low', 'CTG': 'high', 'CTA': 'low', 'CTC': 'high', 'CAC': 'high', 'AAA': 'low', 'CCG': 'high', 'AGT': 'low', 'CCA': 'low', 'CAA': 'low', 'CCC': 'high', 'TAT': 'low', 'GGT': 'low', 'TGT': 'low', 'CGA': 'low', 'CAG': 'high', 'TCT': 'low', 'GAT': 'low', 'CGG': 'high', 'TTT': 'low', 'TGC': 'high', 'GGG': 'high', 'TAG': 'high', 'GGA': 'low', 'TGG': 'high', 'GGC': 'high', 'TAC': 'high', 'TTC': 'high', 'TCG': 'high', 'TTA': 'low', 'TTG': 'high', 'TCC': 'high', 'GAA': 'low', 'TAA': 'high', 'GCA': 'low', 'GTA': 'low', 'GCC': 'high', 'GTC': 'high', 'GCG': 'high', 'GTG': 'high', 'GAG': 'high', 'GTT': 'low', 'GCT': 'low', 'TGA': 'high', 'GAC': 'high', 'CGT': 'low', 'TCA': 'low', 'ATG': 'high', 'CGC': 'high'}

lowAffinityCodons = set(['CTT', 'ACA', 'AAA', 'ATA', 'CCT', 'AGA', 'CAT', 'AAT', 'ATT', 'CTA', 'ACT', 'CAA', 'AGT', 'CCA', 'TAT', 'GGT', 'TGT', 'CGA', 'TCT', 'GAT', 'TTT', 'GGA', 'TTA', 'CGT', 'GAA', 'TCA', 'GCA', 'GTA', 'GTT', 'GCT'])
highAffinityCodons = set(['ACC', 'ATG', 'AAG', 'ACG', 'ATC', 'AAC', 'AGG', 'AGC', 'CTG', 'CTC', 'CAC', 'CCG', 'CAG', 'CCC', 'CGC', 'CGG', 'TGC', 'GGG', 'TAG', 'TGG', 'GGC', 'TAC', 'TTC', 'TCG', 'TTG', 'TCC', 'TAA', 'GCC', 'GTC', 'GCG', 'GTG', 'GAG', 'TGA', 'GAC'])

#</7iyed>

#Empiraclly calculated using genome GRCh37.74 and Ensembl annotations
humanCodonCounts = {'CTT': 588990, 'ACC': 760250, 'ACA': 671093, 'ACG': 248588, 'ATC': 819539, 'AAC': 777291, 'ATA': 326568, 'AGG': 520514, 'CCT': 784233, 'ACT': 581281, 'AGC': 826157, 'AAG': 1373474, 'AGA': 560614, 'CAT': 487348, 'AAT': 745200, 'ATT': 685951, 'CTG': 1579105, 'CTA': 311963, 'CTC': 772503, 'CAC': 618558, 'AAA': 1111269, 'CCG': 285345, 'AGT': 558788, 'CCA': 771391, 'CAA': 572531, 'CCC': 809928, 'TAT': 507376, 'GGT': 459267, 'TGT': 443487, 'CGA': 276584, 'CAG': 1483627, 'TCT': 675336, 'GAT': 982540, 'CGG': 477748, 'TTT': 721642, 'TGC': 495033, 'GGG': 661842, 'TAG': 28685, 'GGA': 731598, 'TGG': 535340, 'GGC': 877641, 'TAC': 588108, 'TTC': 774303, 'TCG': 185384, 'TTA': 348372, 'TTG': 563764, 'TCC': 729893, 'GAA': 1355256, 'TAA': 37503, 'GCA': 718158, 'GTA': 316640, 'GCC': 1120424, 'GTC': 576027, 'GCG': 289438, 'GTG': 1119171, 'GAG': 1685297, 'GTT': 486471, 'GCT': 806491, 'TGA': 82954, 'GAC': 1033108, 'CGT': 200762, 'TCA': 569093, 'ATG': 935789, 'CGC': 404889}

humanCodonCount = 42433513

humanCodonRatios = {'CTT': 0.013880302580651288, 'ACC': 0.017916263496731935, 'ACA': 0.01581516477318293, 'ACG': 0.005858294127097137, 'ATC': 0.019313484603549088, 'AAC': 0.018317856454637634, 'ATA': 0.007695992551924702, 'AGG': 0.012266578070026868, 'CCT': 0.018481453562423644, 'ACT': 0.0136986301369863, 'AGC': 0.01946944623698726, 'AAG': 0.03236767127906662, 'AGA': 0.013211585851965638, 'CAT': 0.011484978865643295, 'AAT': 0.017561591000019253, 'ATT': 0.016165312544356155, 'CTG': 0.0372136287655467, 'CTA': 0.007351807049300867, 'CTC': 0.018205021111497414, 'CAC': 0.014577110313727737, 'AAA': 0.02618847513285077, 'CCG': 0.006724519838835875, 'AGT': 0.01316855382678309, 'CCA': 0.018178815409414725, 'CAA': 0.013492425197037068, 'CCC': 0.01908698909750885, 'TAT': 0.011956964298477951, 'GGT': 0.010823214189218791, 'TGT': 0.010451338308944631, 'CGA': 0.006518055669819277, 'CAG': 0.034963567593378375, 'TCT': 0.015915156494349172, 'GAT': 0.023154811622596506, 'CGG': 0.011258742588670422, 'TTT': 0.017006416602839365, 'TGC': 0.01166608571861585, 'GGG': 0.015597153127529177, 'TAG': 0.0006759987088507143, 'GGA': 0.017241042475083315, 'TGG': 0.012615971720276849, 'GGC': 0.020682732537369696, 'TAC': 0.013859517122704406, 'TTC': 0.01824744041342983, 'TCG': 0.004368811038576985, 'TTA': 0.008209831695999339, 'TTG': 0.013285819630347362, 'TCC': 0.017200861969641778, 'GAA': 0.03193834081095289, 'TAA': 0.0008838061557618385, 'GCA': 0.01692431168732129, 'GTA': 0.007462026535488589, 'GCC': 0.026404224415734798, 'GTC': 0.013574812907901357, 'GCG': 0.006820976618174413, 'GTG': 0.026374695868334068, 'GAG': 0.039716179049328296, 'GTT': 0.011464311239090669, 'GCT': 0.01900599179709679, 'TGA': 0.0019549170958341345, 'GAC': 0.024346511211551115, 'CGT': 0.004731213274752906, 'TCA': 0.013411404330346158, 'ATG': 0.022053064520017467, 'CGC': 0.009541727077840574}

codonTable = {
'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

'CTT' : 'L', 'CTC' : 'L', 'CTA' : 'L', 'CTG' : 'L',
'CCT' : 'P', 'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P',
'CAT' : 'H', 'CAC' : 'H', 'CAA' : 'Q', 'CAG' : 'Q',
'CGT' : 'R', 'CGC' : 'R', 'CGA' : 'R', 'CGG' : 'R',

'ATT' : 'I', 'ATC' : 'I', 'ATA' : 'I', 'ATG' : 'M',
'ACT' : 'T', 'ACC' : 'T', 'ACA' : 'T', 'ACG' : 'T',
'AAT' : 'N', 'AAC' : 'N', 'AAA' : 'K', 'AAG' : 'K',
'AGT' : 'S', 'AGC' : 'S', 'AGA' : 'R', 'AGG' : 'R',

'GTT' : 'V', 'GTC' : 'V', 'GTA' : 'V', 'GTG' : 'V',
'GCT' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A',
'GAT' : 'D', 'GAC' : 'D', 'GAA' : 'E', 'GAG' : 'E',
'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G'
}


AATable = {'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'C': ['TGT', 'TGC'], 'E': ['GAG', 'GAA'], 'D': ['GAT', 'GAC'], 'G': ['GGT', 'GGG', 'GGA', 'GGC'], 'F': ['TTT', 'TTC'], 'I': ['ATC', 'ATA', 'ATT'], 'H': ['CAT', 'CAC'], 'K': ['AAG', 'AAA'], '*': ['TAG', 'TGA', 'TAA'], 'M': ['ATG'], 'L': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'], 'N': ['AAC', 'AAT'], 'Q': ['CAA', 'CAG'], 'P': ['CCT', 'CCG', 'CCA', 'CCC'], 'S': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'], 'R': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'], 'T': ['ACA', 'ACG', 'ACT', 'ACC'], 'W': ['TGG'], 'V': ['GTA', 'GTC', 'GTG', 'GTT'], 'Y': ['TAT', 'TAC']}

AAs = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', '*', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']

toFloat = lambda x: float(x)
toInt = lambda x: int(x)
floatToStr = lambda x:"%f"%(x)
intToStr = lambda x:"%d"%(x)

splitStr = lambda x: x.split(';')
stripSplitStr = lambda x: x.strip().split(';')


def findAll(haystack, needle) :
	"""returns a list of all occurances of needle in haystack"""
	
	h = haystack
	res = []
	f = haystack.find(needle)
	offset = 0
	while (f >= 0) :
		#print h, needle, f, offset
		res.append(f+offset)
		offset += f+len(needle)
		h = h[f+len(needle):]
		f = h.find(needle)

	return res

def reverseComplement(seq):
	'''
	Complements a DNA sequence, returning the reverse complement.
	'''
	return complement(seq)[::-1]

def complement(seq) :
	"""returns the complementary sequence without inversing it"""
	tb = string.maketrans("ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
						  "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn")
	
	return seq.translate(tb)

def translateDNA_6Frames(sequence) :
	"""returns 6 translation of sequence. One for each reading frame"""
	trans = (
				translateDNA(sequence, 'f1'),
				translateDNA(sequence, 'f2'),
				translateDNA(sequence, 'f3'),

				translateDNA(sequence, 'r1'),
				translateDNA(sequence, 'r2'),
				translateDNA(sequence, 'r3'),
			)

	return trans

def translateDNA(sequence, frame = 'f1') :
	"""Translates DNA code, frame : fwd1, fwd2, fwd3, rev1, rev2, rev3"""

	protein = ""

	if frame == 'f1' :
		dna = sequence
	elif frame == 'f2':
		dna = sequence[1:]
	elif frame == 'f3' :
		dna = sequence[2:]
	elif frame == 'r1' :
		dna = reverseComplement(sequence)
	elif frame == 'r2' :
		dna = reverseComplement(sequence)
		dna = dna[1:]
	elif frame == 'r3' :
		dna = reverseComplement(sequence)
		dna = dna[2:]
	else :
		raise ValueError('unknown reading frame: %s, should be one of the following: fwd1, fwd2, fwd3, rev1, rev2, rev3' % frame)

	for i in range(0, len(dna),  3) :
		if (len(dna[i:i+3]) == 3) :
			try :
				protein += codonTable[dna[i:i+3]]
			except KeyError :
				combinaisons = polymorphicCodonCombinaisons(list(dna[i:i+3]))
				translations = set()
				for ci in range(len(combinaisons)):
					translations.add(codonTable[combinaisons[ci]])
				protein += '/'.join(translations)

	return protein

def getSequenceCombinaisons(polymorphipolymorphicDnaSeqSeq, pos = 0) :
	"""Takes a dna sequence with polymorphismes and returns all the possible sequences that it can yield"""

	if type(polymorphipolymorphicDnaSeqSeq) is not types.ListType :
		seq = list(polymorphipolymorphicDnaSeqSeq)
	else :
		seq = polymorphipolymorphicDnaSeqSeq

	if pos >= len(seq) :
		return [''.join(seq)]

	variants = []
	if seq[pos] in polymorphicNucleotides :
		chars = decodePolymorphicNucleotide(seq[pos])
	else :
		chars = seq[pos]#.split('/')

	for c in chars :
		rseq = copy.copy(seq)
		rseq[pos] = c
		variants.extend(getSequenceCombinaisons(rseq, pos + 1))

	return variants

def polymorphicCodonCombinaisons(codon) :
	"""Returns all the possible amino acids encoded by codon"""
	return getSequenceCombinaisons(codon, 0)

def encodePolymorphicNucleotide(polySeq) :
	"""returns a single character encoding all nucletides of polySeq 
	in a single character. PolySeq must have one of the following forms: 
	['A', 'T', 'G'], 'ATG', 'A/T/G'"""
	
	if type(polySeq) is types.StringType :
		if polySeq.find("/") < 0 :
			sseq = list(polySeq)
		else :
			sseq = polySeq.split('/')
			
	else :
		sseq = polySeq
	
	seq = []
	for n in sseq :
		try :
			for n2 in polymorphicNucleotides[n] :
				seq.append(n2)
		except KeyError :
			seq.append(n)
	
	seq = set(seq)
	
	if len(seq) == 4:
		return 'N'
	elif len(seq) == 3 :
		if 'T' not in seq :
			return 'V'
		elif 'G' not in seq :
			return 'H'
		elif 'C' not in seq :
			return 'D'
		elif 'A' not in seq :
			return 'B'
	elif len(seq) == 2 :
		if 'A' in seq and 'G' in seq :
			return 'R'
		elif 'C' in seq and 'T' in seq :
			return 'Y'
		elif 'A' in seq and 'C' in seq :
			return 'M'
		elif 'T' in seq and 'G' in seq :
			return 'K'
		elif 'A' in seq and 'T' in seq :
			return 'W'
		elif 'C' in seq and 'G' in seq :
			return 'S'
	elif polySeq[0] in nucleotides :
		return polySeq[0]
	else :
		raise UnknownNucleotide(polySeq)

def decodePolymorphicNucleotide(nuc) :
	"""the opposite of encodePolymorphicNucleotide, from 'R' to ['A', 'G']"""
	if nuc in polymorphicNucleotides :
		return polymorphicNucleotides[nuc]

	if nuc in nucleotides :
		return nuc

	raise ValueError('nuc: %s, is not a valid nucleotide' % nuc)

def decodePolymorphicNucleotide_str(nuc) :
	"""same as decodePolymorphicNucleotide but returns a string instead 
	of a list, from 'R' to 'A/G"""
	return '/'.join(decodePolymorphicNucleotide(nuc))

def getNucleotideCodon(sequence, x1) :
	"""Returns the entire codon of the nucleotide at pos x1 in sequence, 
	and the position of that nocleotide in the codon in a tuple"""

	if x1 < 0 or x1 >= len(sequence) :
		return None

	p = x1%3
	if p == 0 :
		return (sequence[x1: x1+3], 0)
	elif p ==1 :
		return (sequence[x1-1: x1+2], 1)
	elif p == 2 :
		return (sequence[x1-2: x1+1], 2)

def showDifferences(seq1, seq2) :
	"""Returns a string highligthing differences between seq1 and seq2:
	
	* Matches by '-'
	
	* Differences : 'A|T'
	
	* Exceeded length : '#'
	
	"""
	ret = []
	for i in range(max(len(seq1), len(seq2))) :

		if i >= len(seq1) :
			c1 = '#'
		else :
			c1 = seq1[i]
		if i >= len(seq2) :
			c2 = '#'
		else :
			c2 = seq2[i]

		if c1 != c2 :
			ret.append('%s|%s' % (c1, c2))
		else :
			ret.append('-')

	return ''.join(ret)

def highlightSubsequence(sequence, x1, x2, start=' [', stop = '] ') :
	"""returns a sequence where the subsequence in [x1, x2[ is placed 
	in bewteen 'start' and 'stop'"""

	seq = list(sequence)
	print x1, x2-1, len(seq)
	seq[x1] = start + seq[x1]
	seq[x2-1] = seq[x2-1] + stop
	return ''.join(seq)

# def highlightSubsequence(sequence, x1, x2, start=' [', stop = '] ') :
# 	"""returns a sequence where the subsequence in [x1, x2[ is placed 
# 	in bewteen 'start' and 'stop'"""

# 	seq = list(sequence)
# 	ii = 0
# 	acc = True
# 	for i in range(len(seq)) :
# 		if ii == x1 :
# 			seq[i] = start+seq[i]
# 		if ii == x2-1 :
# 			seq[i] = seq[i] + stop

# 		if i < len(seq) - 1 :
# 			if seq[i+1] == '/':
# 				acc = False
# 			else :
# 				acc = True

# 		if acc :
# 			ii += 1

# 	return ''.join(seq)
