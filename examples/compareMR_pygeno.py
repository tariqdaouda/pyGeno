from pyGeno import Genome
from pyGeno.tools import UsefulFunctions as uf
from pyGeno.Mixed import *
import pyGeno.PositionConverter as pc
from pyGeno.tools import SegmentTree
from filters import *

M=MixedGenome(["human/lightM_transcriptome","human/lightM_Exome"])
R=MixedGenome(["human/lightR_transcriptome","human/lightR_Exome"])
M.loadAllChromosomes()
R.loadAllChromosomes()
nbDiff=0
nbAminoAcids=0
lenNuc = 0

for chroM in M.getChromosomes():
	if chroM.number != 'Y' :
		print 'chro', chroM.number,'...'

		chroR=R.loadChromosome(chroM.number)
		chroM.loadAllGenes(snvsFilter_totalMix_F2) #filterMixed  filterQual
		chroR.loadAllGenes(snvsFilter_totalMix_F2)
		#segment tree : isoformes
		checkSet=set()
		tree = SegmentTree.SegmentTree()
		for geneM in chroM.getGenes():
			geneR=chroR.loadGene(geneM.symbol, snvsFilter_totalMix_F2)
			for transM in geneM.getTranscripts():
				transR=geneR.loadTranscript(transM.id)
				protM=transM.loadProtein() 
				protR=transR.loadProtein()
				
				for e in transR.exons:
					if e.hasCDS():
						tree.insert(e.CDS[0], e.CDS[1], "%s, %s" % (str(e.transcript.id), e.transcript.gene.strand))
						lenNuc += e.getCDSLength()
					
				starFound=False
				i = 0
				while (i < len(protM) and not starFound):
					if "*" in protM[i] or "*" in protR[i]:
						starFound=True
					if protM[i] != protR[i] and not starFound:
						dnaPos = pc.proteinToDNA(i, protM)
						checkSet.add(dnaPos)
					i+=1
				
		nbDiff += len(checkSet)
		nbAminoAcids += tree.getIndexedLength()/3.0
		indLen = tree.getIndexedLength()
		print 'nbAA, nbDiff', nbAminoAcids, nbDiff, nbDiff/float(nbAminoAcids)
		print 'indexLen tree', indLen/3, indLen
		print 'len tree.children', len(tree.children)
		s = 0
		for c in tree.children :
			s+= len(c)
		print "taille tout enfant, lenNuc", s, s/float(indLen), lenNuc
		#print tree
		print "______"
	
	
print "####### final"
print nbDiff
print nbAminoAcids
