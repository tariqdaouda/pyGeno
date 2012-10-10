import string, os, copy, pickle
import numpy as N

class UnknownNucleotide(Exception) :
	def __init__(self, nuc) :
		self.msg =  'Unknown nucleotides %s' % str(nuc)

	def __str__(self) :
		return self.msg
		
nucleotides = ['A', 'T', 'C', 'G']
polymorphicNucleotides = {
			'R' : ['A','G'], 'Y' : ['C','T'], 'M': ['A','C'],
			'K' : ['T','G'], 'W' : ['A','T'], 'S' : ['C','G'], 'B': ['C','G','T'],
			'D' : ['A','G','T'], 'H' : ['A','C','T'], 'V' : ['A','C','G'], 'N': ['A','C','G','T']
			}

#from Molecular Systems Biology 8; Article number 572; doi:10.1038/msb.2012.3
lowAffinityCodons = ['GCA', 'GCT', 'AGA', 'CGA', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGA', 'GGT', 'CAT', 'ATA', 'ATT', 'CTA', 'CTT', 'TTA', 'AAA', 'TTT', 'CCA', 'CCT', 'AGT', 'TCA', 'TCT', 'ACA', 'ACT', 'TAT', 'GTA', 'GTT']

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

toFloat = lambda x: float(x)
toInt = lambda x: int(x)
floatToStr = lambda x:"%f"%(x)
intToStr = lambda x:"%d"%(x)

splitStr = lambda x: x.split(';')
stripSplitStr = lambda x: x.strip().split(';')

"""
def strToList(stri) :
	l = []
	for c in seqtri :
		l.append(c)
	
	return l
"""	
def saveResults(directoryName, results, errors = '', params = '', fileName = ''):
	
	if not os.path.exists("results/" + directoryName):
		os.makedirs("results/" + directoryName)
	
	resPath = "results/%s/%s"%(directoryName, fileName)
	resFile = open(resPath, 'w')
	print "Saving results :\n\t%s..."%resPath
	resFile.write(results)
	resFile.close()

	if errors != '' :
		errPath = "results/%s/%s.err.txt"%(directoryName, fileName)
		errFile = open(errPath, 'w')

		print "Saving error log :\n\t%s..." %errPath
		errFile.write(errors)
		errFile.close()
	
	if params != '' :
		paramPath = "results/%s/%s.params.txt"%(directoryName, fileName)
		paramFile = open(paramPath, 'w')

		print "Saving params :\n\t%s..." %paramPath
		paramFile.write(params)
		paramFile.close()
		
	return "results/%s/"%(directoryName)

def saveResults_pkl(filename, obj):
	if filename.find('.pygeno-pklresults') > 1:
		f = open(filename, 'w')
	else :
		f = open(filename+'.pygeno-pklresults', 'w')

	pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
	f.close()

def unescapeHtml(html) :
	h = html
	h = h.replace('&lt;', '<')
	h = h.replace('&gt;', '>')
	h = h.replace('&amp;', '&')
	h = h.replace('&quot;', '"')
	h = h.replace('&apos;', "'")

	return h

def removeHtmlFromEfetch(html) :
	pos = html.find('<pre>')
	fpos = html.find('</pre>')
	return html[pos+5:fpos].strip()

def toXML(template, data, addHeader = False):
	"""
	@param template
		<snpDistribution>
			<geneSymbol>#GENESYMB#</geneSymbol>
			<proteinId>#PROTID#</proteinId>
			<distribution>#DIST#</distribution>
			#SNPLIST#
		</snpDistribution>
	@param data : data['#GENESYMB#'] = xxxx;...
	@param addHeader : appends <?xml version="1.0" encoding="UTF-8"?> at begining of the file
	@returns : data formated in an xml fashion"""
	t = template
	for k in data.keys():
		t = t.replace(k, data[k])
	
	if addHeader :
		return '<?xml version="1.0" encoding="UTF-8"?>\n\t'+t

	return t


"""returns a list of all the ocurances"""
def findAll(haystack, needle) :
	
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
	"returns the complementary sequence without inversing it"
	tb = string.maketrans("ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
						  "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn")
	return seq.translate(tb)
  
def translateDNA_6Frames(sequence) :
	trans = [
				translateDNA(sequence, 1),
				translateDNA(sequence, 2),
				translateDNA(sequence, 3),
				
				translateDNA(sequence, -1),
				translateDNA(sequence, -2),
				translateDNA(sequence, -3),
			]

	return trans

"""Translates DNA code, frame : 1, 2, 3, -1, -2, -3"""
def translateDNA(sequence, frame = 1) :

	protein = ""
	
	if frame == 1 :
		dna = sequence
	elif frame == 2:
		dna = sequence[1:]
	elif frame == 3 :
		dna = sequence[2:]
	elif frame == -1 :
		dna = complement(sequence)
	elif frame == -2 :
		dna = complement(sequence)
		dna = dna[1:]
	elif frame == -3 :
		dna = complement(sequence)
		dna = dna[2:]

	for i in range(0, len(dna),  3) :
		if (len(dna[i:i+3]) == 3) :
			try :
				protein += codonTable[dna[i:i+3]]
				#print i, dna[i:i+3], codonTable[dna[i:i+3]]
			except KeyError :
				#print dna[i:i+3]
				combinaisons = polymorphicCondonCombinaisons(list(dna[i:i+3]))
				translations = set()
				for ci in range(len(combinaisons)):
					translations.add(codonTable[combinaisons[ci]])
				protein += '/'.join(translations)
				
	return protein
	
def polymorphicCondonCombinaisons(dnaSeq, startId = 0) :
	if type(dnaSeq).__name__ != 'list' :
		dna = list(dnaSeq)
	else :
		dna = dnaSeq
	
	if startId >= len(dna) :
		return [''.join(dna)]
		
	rDna = copy.copy(dna)
	res = []
	try :
		chars = polymorphicNucleotides[dna[startId]]
		for c in chars :
			rDna[startId] = c
			res.extend(polymorphicCondonCombinaisons(rDna, startId +1))
		
	except KeyError:
		res.extend(polymorphicCondonCombinaisons(rDna, startId +1))

	return res
	
def getPolymorphicNucleotide(strSeq) :
	"""Seq is a string like : ATG"""
	
	seq = list(strSeq)
	for i in range(len(seq)) :
		if seq[i] in polymorphicNucleotides :
			seq[i] = ''.join(polymorphicNucleotides[seq[i]])
	
	seq = set(strSeq)
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
	elif strSeq[0] in nucleotides :
		return strSeq[0]
	else :
		raise UnknownNucleotide(strSeq)

def decodePolymorphicNucleotide(nuc) :
	if nuc in polymorphicNucleotides :
		return '/'.join(polymorphicNucleotides[nuc])
	
	return nuc
	
def getCodon(sequence, x1) :
	"Returns the entire codon of the nucleotide at pos x1 in the cdna, and the position of that nocleotide in the codon"
	
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
	"returns a string representig matching chars (-) and differences (A/T) between the two strings or length exceeded (#)"
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
			ret.append('%s/%s' % (c1, c2))
		else :
			ret.append('-')
		
	return ''.join(ret)
	
def highLightSubsequence(sequence, x1, x2, start=' [', stop = '] ') :
	"returns a sequence where the subsequence in [x1, x2[ is placed bewteen start and stop,"
	
	seq = list(sequence)
	ii = 0
	acc = True
	for i in range(len(seq)) :
		if ii == x1 :
			seq[i] = '['+seq[i]
		if ii == x2-1 :
			seq[i] = seq[i] + ']'
		
		if i < len(seq-1) :
			if seq[i+1] == '/':
				acc = False
			else :
				acc = True
				
		if acc :
			ii += 1
			
	print ''.join(seq)
