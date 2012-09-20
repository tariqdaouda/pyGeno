from tools import UsefulFunctions as uf
import copy

#class PositionConverter :

#	def __init__(self) :
#		pass

def cDNAToDNA(position, transcript, iWantCDSNumber = False) :
	"""Returns -1 if pos is outside the transcript
	if iWantCDSNumber = returns a tuple (position, cds), if false returns only the postion"""
	#if transcript.gene.strand == '-' :
		#pos = len(transcript) - position
		#CDSs = transcript.CDSs
	#else :
	pos = position
		#CDSs = transcript.CDSs[::-1]
	
	#print "----", position, len(transcript), transcript.getNumberOfCDSs()
	
	poffset = pos
	prevLen = 0
	#i = 0
	
	#for cds in transcript.CDSs :
	for i in range(len(transcript.exons)) :
		cds = transcript.exons[i].CDS
		
		if cds != None :
			if transcript.gene.strand == '+' :
				poffset = (pos - prevLen) + cds[0]
			else :
				poffset = cds[1] - (pos - prevLen)

			#print poffset, cds.x2
			#print "-8-8-8-", pos, poffset, cds.number, (cds.x1, cds.x2), (cds.x1 <= poffset),  (poffset <= cds.x2)
			if (cds[0] <= poffset) and (poffset <= cds[1]) :
				#print transcript.gene.chromosome.getSequence(poffset, poffset+1)
				#sa marche mais je ne sais pas pourquoi
				#if transcript.gene.strand == '-' :
				#	poffset += 1
				
				if iWantCDSNumber :
					return (poffset, i)
				return poffset
			
			prevLen += transcript.exons[i].getCDSLength() + 1 #+1 premiere nuc du nouvo cds
		#if transcript.gene.strand == '+' :
			
		#else :
		#	prevLen -= len(cds)
		#i+=1	
		
	return -1


def cDNAToDNA_range(x1, x2, transcript, iWantCDSNumber = False) :

	xx1 = cDNAToDNA(x1, transcript, True)
	xx2 = cDNAToDNA(x2, transcript, True)

	#print '\t', xx1, xx2
	if xx1 == -1 or xx2 == -1 :
		return None
	
	if xx2[0] < xx1[0] :
		#print "inv", xx1[1].number,  xx2[1].number
		tmp = xx1
		xx1 = xx2
		xx2 = tmp

	if xx1[1] == xx2[1] :			
		if iWantCDSNumber :
			ret = [(xx1[0], xx2[0], xx1[1])]
		else:
			ret = [(xx1[0], xx2[0])]
	else :
		ret = []
		
		cds1 = transcript.exons[xx1[1]].CDS
		cds2 = transcript.exons[xx2[1]].CDS
	
		if cds1 != None and cds2 != None :
			if iWantCDSNumber :
				ret.append((xx1[0], cds1[1], xx1[1]))
			else :
				ret.append((xx1[0], cds1[1]))
				
			for i in range(xx1[1]+1, xx2[1]):
				#print range(xx1[1].position+1, xx2[1].position) 
				cdsTemp = transcript.exons[i].CDS
				if cdsTemp != None :
					if iWantCDSNumber :
						ret.append((cdsTemp[0], cdsTemp[1], i))
					else:
						ret.append((cdsTemp[0], cdsTemp[1]))
			
			#print 'bbb', cds2.x1, xx2[0]
			if iWantCDSNumber :
				ret.append((cds2[0], xx2[0], xx2[1]))
			else :
				ret.append((cds2[0], xx2[0]))
		
	#for a, b in ret :
	#	print 'aaa', a, b
	return ret
	
def DNAToCDNA(pos, transcript) :
	"""Returns -1 if the position is outside the geneTODO a tester"""
	
	prevLen = 0
	for cds in transcript.CDSs :
		pres = pos - cds.x1 + prevLen
		#print "===>", pres, pos, (cds.x1, cds.x2), cds.number
		if (cds.x1 <= pos) and (pos <= cds.x2) :
			if transcript.gene.strand == '-' :
				return len(transcript) - 1 - pres
			else :
				return pres
				
		prevLen += len(cds)
		
	return -1

def proteinToCDNA(pos) :
	return int(pos*3)
		
def cDNAtoProtein(pos) :
	return int(pos)/3

def proteinToDNA(pos, protein) :
	return cDNAToDNA(proteinToCDNA(pos), protein.transcript) 

def proteinToDNA_range(x1, x2, protein, iWantCDSNumbers = False) :
	xx1 = proteinToCDNA(x1)
	xx2 = proteinToCDNA(x2)+2 #x2 est la pos de la premiere nuc d'ou +2
	#print xx1, xx2
	return cDNAToDNA_range(xx1, xx2, protein.transcript, iWantCDSNumbers) 
	
def DNAtoProtein(pos, protein):
	return cDNAtoProtein(DNAToCDNA(pos, protein.transcript))
	
