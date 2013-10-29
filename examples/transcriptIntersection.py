from pyGeno.tools.SegmentTree import SegmentTree
from pyGeno.Genome import Genome
import pyGeno.PositionConverter as PC
from pyGeno.tools import UsefulFunctions as uf

print "Example: print the cdna of all transcripts that intersect given positions. print the corresponding codon and the amino acid"

positions = [77170087, 77170095]

genome = Genome('human/reference')
chro = genome.loadChromosome('18')

index = SegmentTree()

print "building index for chromosome %s..." % chro.number
for gene in chro.getGenes():
	print "\t", gene.symbol
	for trans in gene.getTranscripts():
		for exon in trans.getExons() :
			
			#To only index the coding sequences
			if exon.hasCDS() :
				index.insert(exon.CDS[0], exon.CDS[1], name = str(exon), referedObject = [trans])
				#print exon.CDS
				
			#the whole exon
			#index.insert(exon.x1, exon.x2, name = str(exon), referedObject = [trans])
			
			#You may want to pickle the index, so you don't have to regenerate it at each time
			#This can result in a pretty large object due to the chain references. One way to avoid that
			#problem is to use trans.pluck() as the referedObject. By plucking the object you remove
			#all chain references, effectively saving the transcript alone
			
print "intersections..."
for pos in positions :
	#get the list of segments intersected by pos
	interections = index.intersect(pos)
	for inter in interections :
		for trans in inter.referedObject :
			print "-----"
			print trans.CDNA #trans.sequence for the whole sequence
			cdnaPos = PC.DNAToCDNA(pos, trans)
			codon = trans.getCodon(cdnaPos)

			print "dna position %s, cdna position %s, codon %s, aa %s" %(pos, cdnaPos, codon, uf.translateDNA(codon[0]))
