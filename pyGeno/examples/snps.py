from pyGeno.Genome import Genome
from pyGeno.Gene import Gene
from pyGeno.Transcript import Transcript

ref = Genome(name = 'GRCh37.75')
dummy = Genome(name = 'GRCh37.75', SNPs = 'dummySRY')
geneRef = ref.get(Gene, name = 'SRY')[0]
for trans in geneRef.get(Transcript) :
	print trans.name
	print trans.protein.sequence
	for e in trans.exons :
		print e.CDS, e.strand,  e.CDS_start, e.CDS_end
print '-------'
geneDummy = dummy.get(Gene, name = 'SRY')[0]
for trans in geneDummy.get(Transcript) :
	print trans.name
	print trans.protein.sequence
	for e in trans.exons :
		print e.CDS, e.strand,  e.CDS_start, e.CDS_end
