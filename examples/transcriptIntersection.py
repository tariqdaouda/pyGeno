from pyGeno.tools.SegmentTree import SegmentTree
from pyGeno.Genome import Genome

print "Example: print the sequences of all transcripts that intersect given posotions"

positions = [2982666, 2982870]

genome = Genome('human/reference')
chro = genome.loadChromosome('18')

index = SegmentTree()

print "building index for chromosome %s..." % chro.number
for gene in chro.getGenes():
	print "\t", gene.symbol
	for trans in gene.getTranscripts():
		for exon in trans.getExons() :
			
			index.insert(exon.x1, exon.x2, name = str(exon), referedObject = [trans])
			#To only index the coding sequences
			#if exon.hasCDS :
			#	index.insert(exon.CDS[0], exon.CDS[1], name = str(exon), referedObject = [trans])
			
			#You may want to pickle the index, so you don't have to regenerate it at each time
			#This can result in a pretty large objects due to the chain references. One way to avoir that
			#problem is to use trans.pluck() as the referedObject. By plucking the object you remove
			#all chain references effectively saving the transcript alone
			
print "intersections..."
for pos in positions :
	#retourne une liste de tout les segment intersectes par pos
	interections = index.intersect(pos)
	for inter in interections :
		for trans in inter.referedObject :
			print trans.sequence # trans.CDNA pour la sequence codante
