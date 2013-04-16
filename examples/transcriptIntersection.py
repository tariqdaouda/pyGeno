from pyGeno.tools.SegmentTree import SegmentTree
from pyGeno.Genome import Genome
import pyGeno.PositionConverter as PC

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
			#print exon.x1,  exon.x2
			#ou si on veux juste les sequences codantes
			#if exon.hasCDS :
			#	index.insert(exon.CDS[0], exon.CDS[1], name = str(exon), referedObject = [trans])

print "intersections..."
for pos in positions :
	#retourne une liste de tout les segment intersectes par pos
	interections = index.intersect(pos)
	for inter in interections :
		for trans in inter.referedObject :
			print trans.sequence # trans.CDNA pour la sequence codante
