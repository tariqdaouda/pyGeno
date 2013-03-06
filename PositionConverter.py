from tools import UsefulFunctions as uf
import copy

def CDNAToDNA(position, transcript, iWantCDSNumber = False) :
	"""Returns -1 if pos is outside the transcript
	if iWantCDSNumber = returns a tuple (position, cds), if false returns only the postion"""

	pos = position		
	poffset = 0

	for i in range(len(transcript.codingExons)) :
		cds = transcript.codingExons[i].CDS

		if pos < transcript.codingExons[i].getCDSLength() or (i == len(transcript.codingExons)-1 and pos == transcript.codingExons[i].getCDSLength()) :
			
			if transcript.gene.strand == '+' :
				poffset = cds[0] + pos
			else :
				poffset = cds[1] - pos-1
		
			if iWantCDSNumber :
				#print '---', position, poffset, transcript.codingExons[i].CDS
				return (poffset,i)
			return poffset
		
		pos -= transcript.codingExons[i].getCDSLength()

	return -1


def CDNAToDNA_range(x1, x2, transcript, iWantCDSNumber = False) :
	
	xx1 = CDNAToDNA(x1, transcript, True)
	xx2 = CDNAToDNA(x2, transcript, True)
	
	#print 'xx a', xx1, xx2
	if xx2[1] < xx1[1] :
		tmp = xx2
		xx2 = xx1
		xx1 = xx2
	
	if transcript.gene.strand == '-' :
		xx1 = (xx1[0]+1, xx1[1])
		xx2 = (xx2[0]+1, xx2[1])
	
	#print 'xx b', xx1, xx2
	
	#print xx1, xx2
	if xx1[1] == xx2[1] :
		if iWantCDSNumber :
			return [(xx1[0], xx2[0], xx1[1])]
		else :
			return [(xx1[0], xx2[0])]
	
	if transcript.gene.strand == '+' :
		if iWantCDSNumber :
			ret = [(xx1[0], transcript.codingExons[xx1[1]].CDS[1], xx1[1])]
		else :
			ret = [(xx1[0], transcript.codingExons[xx1[1]].CDS[1])]
	else :
		if iWantCDSNumber :
			ret = [(xx1[0], transcript.codingExons[xx1[1]].CDS[0], xx1[1])]
		else :
			ret = [(xx1[0], transcript.codingExons[xx1[1]].CDS[0])]
			
	for i in range(xx1[1]+1, xx2[1]) :
		cds = transcript.codingExons[i].CDS
		if iWantCDSNumber :
			ret.append((cds[0], cds[1], i))
		else :
			ret.append((cds[0], cds[1]))
	
	cds = transcript.codingExons[xx2[1]].CDS
	if transcript.gene.strand == '+' :
		if iWantCDSNumber :
			ret.append((cds[0], xx2[0], i))
		else :
			ret.append((cds[0], xx2[0]))
	else :
		if iWantCDSNumber :
			ret.append((cds[1], xx2[0], i))
		else :
			ret.append((cds[1], xx2[0]))
	
	if transcript.gene.strand == '-' :
		return ret[::-1]
	return ret
	
def CDNAToDNA_range_bck(x1, x2, transcript, iWantCDSNumber = False) :

	xx1 = CDNAToDNA(x1, transcript, True)
	xx2 = CDNAToDNA(x2, transcript, True)
	if transcript.gene.strand == '-' :
		xx1 = (xx1[0]+1, xx1[1])
		xx2 = (xx2[0]+1, xx2[1])

	if xx1 == -1 or xx2 == -1 :
		print '\t cdna=>dna err', x1, xx1, x2, xx2, transcript, transcript.gene.strand
		return None
	
	
	if xx1[1] == xx2[1] :			
		if iWantCDSNumber :
			ret = [(xx1[0], xx2[0], xx1[1])]
		else:
			ret = [(xx1[0], xx2[0])]
	else :
		ret = []
		
		cds1 = transcript.codingExons[xx1[1]].CDS
		cds2 = transcript.codingExons[xx2[1]].CDS

		if transcript.gene.strand == '+' :
			if iWantCDSNumber :
				#ret.append((xx1[0], cds1[1]+1, xx1[1]))
				ret.append((xx1[0], cds1[1], xx1[1]))
			#ret.append((xx1[0], cds1[1]+1))
			ret.append((xx1[0], cds1[1]))
		else :			
			if iWantCDSNumber :
				ret.append((cds1[0], xx1[0], xx1[1]))
			ret.append((cds1[0], xx1[0]))
			
			
		#for i in range(xx1[1]+1, xx2[1]):
		for i in range(xx1[1], xx2[1]):
			cdsTemp = transcript.codingExons[i].CDS
		
			if iWantCDSNumber :
				#ret.append((cdsTemp[0], cdsTemp[1]+1, i))
				ret.append((cdsTemp[0], cdsTemp[1], i))
			#ret.append((cdsTemp[0], cdsTemp[1]+1))
			ret.append((cdsTemp[0], cdsTemp[1]))
	
		if transcript.gene.strand == '+' :
			if iWantCDSNumber :
				ret.append((cds2[0], xx2[0], xx2[1]))
			ret.append((cds2[0], xx2[0]))
		else :			
			if iWantCDSNumber :
				#ret.append((xx2[0], cds2[1]+1, xx2[1]))
				ret.append((xx2[0], cds2[1], xx2[1]))
			#ret.append((xx2[0], cds2[1]+1))
			ret.append((xx2[0], cds2[1]))

	if transcript.gene.strand == '-' :
		return ret[::-1]
	return ret
	
def DNAToCDNA(pos, transcript) :
	"return None if the position is out of range"
	cdsLen = 0
	for e in transcript.exons :
		if e.hasCDS() :
			print pos, e.CDS
			if e.CDS[0] <= pos and pos < e.CDS[1]:
				if transcript.gene.strand == '+'  :
					resPos = pos - e.CDS[0] + cdsLen
					return resPos
				resPos = (e.CDS[1] -pos) + cdsLen -1
				return resPos
			else :
				cdsLen += e.getCDSLength()
	return None
	
def proteinToCDNA(pos) :
	return int(pos*3)

def CDNAToProtein(pos) :
	return int(pos)/3

def proteinToDNA(pos, protein) :
	return CDNAToDNA(proteinToCDNA(pos), protein.transcript) 

def proteinToDNA_range(x1, x2, protein, iWantCDSNumbers = False) :
	xx1 = proteinToCDNA(x1)
	xx2 = proteinToCDNA(x2)
	
	return CDNAToDNA_range(xx1, xx2, protein.transcript, iWantCDSNumbers) 
	
def DNAToProtein(pos, protein):
	return CDNAToProtein(DNAToCDNA(pos, protein.transcript))

def chromosomeToGenome(pos, chromosome):
	if pos > len(chromosome) or pos < 0:
		raise IndexError('%d not in chromosome %s of len %d' % (pos, chromosome.number, len(chromosome)))

	return pos + chromosome.x1
	
def genomeToChromosome(pos, chromosome):
	if pos > len(chromosome.genome) or pos < 0:
		raise IndexError('%d not in genome %s of len %d' % (pos, chromosome.genome.name, len(chromosome.genome)))

	p = pos - chromosome.x1
	if p < 0 :
		raise IndexError('%d not in chromosome %s of len %d' % (pos, chromosome.number, len(chromosome.genome)))
		
	return pos - chromosome.x1
