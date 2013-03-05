from pyGeno.Genome import Genome
from pyGeno.tools import UsefulFunctions as uf
from pyGeno.tools.SegmentTree import SegmentTree

def countAAs(seq) :
	counts = {}
	for aa in uf.AAs :
		counts[aa] = 0
	
	for c in seq :
		if c != '/' :
			counts[c] += 1
	
	return counts

resCounts = {}
for aa in uf.AAs :
	resCounts[aa] = 0

B6 = Genome('mouse/B6') 
print 'Loading data for genome %s...' % B6.name
B6.loadAllChromosomes(False)
print 'Fight!'
for chro in B6.chromosomes.keys() :
	chroB6 = B6.loadChromosome(chro)
	chroB6.loadAllGenes()
	for gene in chroB6.genes.keys() :
		geneB6 = chroB6.loadGene(gene)
		for trans in geneB6.transcripts.keys() :
			protB6 = geneB6.loadTranscript(trans).loadProtein()
			counts = countAAs(protB6.sequence)
			for aa in counts :
				resCounts[aa] += counts[aa]

total = 0
for aa in resCounts:
	total += resCounts[aa]

print 'Usage for genome %s' % B6.name

for aa in resCounts:
	print 'aa', aa, 'count', resCounts[aa], 'sur', total, 'ratio', float(resCounts[aa])/total
