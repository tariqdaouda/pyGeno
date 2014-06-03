from pyGeno.Genome import Genome
from pyGeno.Gene import Gene
#from pyGeno.Chromosome import Chromosome
from pyGeno.Transcript import Transcript
#from pyGeno.Protein import Protein

import pyGeno.configuration as conf

#conf.db.enableDebug(True)
print Transcript.help()

refGenome = Genome(name = "GRCh37.75")
trans = refGenome.get(Transcript, name = 'TPST2-001')[0]
print trans.exons[0].sequence
