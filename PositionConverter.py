from tools import UsefulFunctions as uf
import copy

def CDNAToDNA(position, transcript, iWantCDSNumber = False) :
	"""Returns -1 if pos is outside the transcript
	if iWantCDSNumber = returns a tuple (position, cds), if false returns only the postion"""

	pos = position		
	poffset = 0

	for i in range(len(transcript.codingExons)) :
		cds = transcript.codingExons[i].CDS
		 
		#if cds != None :
		#print '===> pos', pos, cds, transcript.codingExons[i].getCDSLength(), transcript.codingExons[i].number, i
		if pos < transcript.codingExons[i].getCDSLength() or (i == len(transcript.codingExons)-1 and pos == transcript.codingExons[i].getCDSLength()) :
			
			if transcript.gene.strand == '+' :
				poffset = cds[0] + pos
			else :
				poffset = cds[1] - pos
			
			#print "poffset", poffset
			#print transcript.gene.chromosome.getSequence(cds[0], cds[1])
			#print uf.complement(transcript.gene.chromosome.getSequence(cds[0], cds[1]))
			if iWantCDSNumber :
				return (poffset,i)
			return poffset
		
		pos -= transcript.codingExons[i].getCDSLength()

	
	return -1


def CDNAToDNA_range(x1, x2, transcript, iWantCDSNumber = False) :

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
		#print "CDS1 !!! aaa", xx1, cds1

		if transcript.gene.strand == '+' :
			if iWantCDSNumber :
				ret.append((xx1[0], cds1[1]+1, xx1[1]))
			ret.append((xx1[0], cds1[1]+1))
		else :			
			if iWantCDSNumber :
				ret.append((cds1[0], xx1[0], xx1[1]))
			ret.append((cds1[0], xx1[0]))
			
			
		for i in range(xx1[1]+1, xx2[1]):
			cdsTemp = transcript.codingExons[i].CDS
			#if transcript.gene.strand == '+' :
			#	if iWantCDSNumber :
			#		ret.append((cdsTemp[0], cdsTemp[1]+1, i))
			#	ret.append((cdsTemp[0], cdsTemp[1]+1, i))
			#else :			
			if iWantCDSNumber :
				ret.append((cdsTemp[0], cdsTemp[1]+1, i))
			ret.append((cdsTemp[0], cdsTemp[1]+1))
	
		if transcript.gene.strand == '+' :
			#print "CDS2 !!! aaa", xx2, cds2
			if iWantCDSNumber :
				ret.append((cds2[0], xx2[0], xx2[1]))
			ret.append((cds2[0], xx2[0]))
		else :			
			if iWantCDSNumber :
				ret.append((xx2[0], cds2[1]+1, xx2[1]))
			ret.append((xx2[0], cds2[1]+1))

	if transcript.gene.strand == '-' :
		return ret[::-1]
	#print '-----', ret, transcript.gene.strand, iWantCDSNumber
	return ret
	
def DNAToCDNA(pos, transcript) :
	"return None if the position is out of range"
	cdsLen = 0
	for e in transcript.exons :
		if e.hasCDS() :
			#print e.CDS, pos, e.CDS[0] <= pos and pos <= e.CDS[1], transcript.gene.strand
			if e.CDS[0] <= pos and pos <= e.CDS[1]:
				if transcript.gene.strand == '+'  :
					resPos = pos - e.CDS[0] + cdsLen
					#print 'respos', resPos, pos, '-', e.CDS[0], '+', cdsLen
					return resPos
				resPos = (e.CDS[1] -pos) + cdsLen #+1
				return resPos
			else :
				cdsLen += e.getCDSLength()
	#print 'NF'
	return None

def proteinToCDNA(pos) :
	return int(pos*3)
		
def CDNAtoProtein(pos) :
	return int(pos)/3

def proteinToDNA(pos, protein) :
	return CDNAToDNA(proteinToCDNA(pos), protein.transcript) 

def proteinToDNA_range(x1, x2, protein, iWantCDSNumbers = False) :
	xx1 = proteinToCDNA(x1)
	xx2 = proteinToCDNA(x2)#+2 #x2 est la pos de la premiere nuc d'ou +2
	#print '==>', x1, xx1,'|',x2, xx2, '|', len(protein.transcript.CDNA), len(protein)*3, protein.transcript.id, protein.transcript.gene.symbol, protein.transcript.gene.chromosome.number
	#print "====>", protein
	return CDNAToDNA_range(xx1, xx2, protein.transcript, iWantCDSNumbers) 
	
def DNAtoProtein(pos, protein):
	return CDNAtoProtein(DNAToCDNA(pos, protein.transcript))
	
