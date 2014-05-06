from pyGeno.Genome import Genome
from pyGeno.Gene import Gene
#from pyGeno.Chromosome import Chromosome
from pyGeno.Transcript import Transcript
#from pyGeno.Protein import Protein

import pyGeno.configuration as conf

#conf.db.enableDebug(True)
print Transcript.help()

refGenome = Genome(name = "GRCh37.74", specie = 'human')
gene = refGenome.get(Gene, name = 'TPST2')[0]
for t in gene.get(Transcript) :
	print t.sequence
