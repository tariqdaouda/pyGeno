from pyGeno.Genome import Genome
from pyGeno.Gene import Gene
#from pyGeno.Chromosome import Chromosome
from pyGeno.Transcript import Transcript
from pyGeno.Exon import *
#from pyGeno.Protein import Protein

import pyGeno.configuration as conf

if False :
	#conf.db.enableDebug(True)
	print Transcript.help()

	refGenome = Genome(name = "GRCh37.75")
	#if False :
	trans = refGenome.get(Transcript, name = 'TPST2-001')[0]
	print trans.sequence
	print trans.exons[0].sequence


	#from pyGeno.Genome import *
	#from pyGeno.Chromosome import *

	e = Exon(raba_id=306123)
	e1 = e.previousExon()
	print e.sequence == e1.sequence
	print e.sequence
	print e1.sequence

	from pyGeno.Genome import *
	from pyGeno.Chromosome import *
	from pyGeno.Exon import *
e = Exon(raba_id=306123)

print e.number
print e.previousExon().number
print e.nextExon().number
