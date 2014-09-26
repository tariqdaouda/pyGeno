pyGeno: A python package for Personalized Proteogenomics
========================================================

Importing a data wrap:
----------------------
from pyGeno.importation.Genomes import *
importGenome('GRCh37.75.tar.gz')

Instanciating a genome:
-----------------------
from pyGeno.Genome import Genome
ref = Genome(name = 'GRCh37.75')

Printing all the proteins of a gene:
-----------------------------------
from pyGeno.Genome import Genome
from pyGeno.Gene import Gene
from pyGeno.Protein import Protein
ref = Genome(name = 'GRCh37.75')
gene = ref.get(Gene, name = 'TPST2')[0]
for prot in gene.get(Protein) :
	print prot.sequence

Progress Bar:
-------------

from pyGeno.tools.ProgressBar import ProgressBar
pg = ProgressBar(nbEpochs = 155)
for i in range(155) :
	p.update(label = '%d' %i) # or simply p.update() 
p.close()

