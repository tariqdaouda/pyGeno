pyGeno: A python package for Personalized Proteogenomics
========================================================
.. code:: python
	
	pip install pyGeno #for the latest stable version

For the latest developments:

.. code:: shell

	git clone https://github.com/tariqdaouda/pyGeno.git
	cd pyGeno
	python setup.py develop

To run tests:

.. code:: shell

	python setup.py test
	
The full documentation is available here_:

.. _here: http://bioinfo.iric.ca/~daoudat/pyGeno/

Bootstraping:
-------------

pyGeno comes with a few datawraps, to get the list you can use:

.. code:: python
	
	import pyGeno.bootstrap as B
	B.listDataWraps()

Importing whole genomes is a demanding process that takes time and requires (according to tests) at least 2GB
of memory. Depending on your configuration more might be required.

That being said importating a datawrap is a one time operation and once the importation is complete the datawrap
can be discarded without consequences.

The bootstrap module also has some handy functions for importing built-in packages.

Some of them just for playing around with pyGeno (**Fast importation** and **Small memory requirements**):

.. code:: python
	
	import pyGeno.bootstrap as B
	
	#Imports only the first and Y chromosomes from the human reference genome GRCh37.75
	#Fast, and does not require much memory. Sequences of both chromosomes will be downloaded. 
	B.importHumanReference_1YOnly()

	#Imports only the Y chromosome from the human reference genome GRCh37.75
	#Very fast, requires even less memory. No download required.
	B.importHumanReference_YOnly()
	
	#A dummy datawrap for humans that mimics a casava's snps.txt with one SNP at the begining of the gene SRY
	B.importDummySRY()

And for more **Serious Work**, the whole reference genome.

.. code:: python

	#Downloads the whole genome (205MB, sequences + annotations), may take an hour or more.
	B.importHumanReference()
	
Importing a data wrap:
----------------------

.. code:: python

  from pyGeno.importation.Genomes import *
  importGenome('GRCh37.75.tar.gz')

To import a patient's specific polymorphisms

.. code:: python

  from pyGeno.importation.SNPs import *
  importSNPs('patient1.tar.gz')

You can easily make your own datawraps with any tar.gz compressor.
For more details on how data wraps are made you can have a look inside the folder bootstrap_data, 
of check the documentation Importation_.

.. _Importation: http://bioinfo.iric.ca/~daoudat/pyGeno/importation.html

Instanciating a genome:
-----------------------
.. code:: python
	
	from pyGeno.Genome import Genome
	#the name of the genome is defined inside the package's manifest.ini file
	ref = Genome(name = 'GRCh37.75')

Printing all the proteins of a gene:
-----------------------------------
.. code:: python

  from pyGeno.Genome import Genome
  from pyGeno.Gene import Gene
  from pyGeno.Protein import Protein
  ref = Genome(name = 'GRCh37.75')
  #get returns a list of elements
  gene = ref.get(Gene, name = 'TPST2')[0]
  for prot in gene.get(Protein) :
  	print prot.sequence

Making queries, get() Vs iterGet():
-----------------------------------
iterGet is a faster version of get that returns an iterator instead of a list.

Making queries, syntax:
----------------------
pyGeno's get function uses the expressivity of rabaDB.

These are all possible query formats:

.. code:: python

  ref.get(Gene, name = "SRY")
  ref.get(Gene, { "name like" : "HLA"})
  chr12.get(Exon, { "start >=" : 12000, "end <" : 12300 })
  ref.get(Transcript, { "gene.name" : 'SRY' })


Creating indexes to speed up queries:
------------------------------------
.. code:: python

  from pyGeno.Gene import Gene
  #creating an index on gene names if it does not already exist
  Gene.ensureGlobalIndex('name')
  #removing the index
  Gene.dropIndex('name')
  
Creating a Personalized Genome:
-------------------------------
Personalized Genomes are a powerful feature that allow to work on the specific genomes and proteomes of your patients.
You can even mix several SNPs together.

.. code:: python
  
  from pyGeno.Genome import Genome
  #the name of the snp set is defined inside the package's manifest.ini file
  dummy = Genome(name = 'GRCh37.75', SNPs = 'dummySRY')
  #you can also define a filter (ex: a quality filter) for the SNPs
  dummy = Genome(name = 'GRCh37.75', SNPs = 'dummySRY', SNPFilter = myFilter)
  #and even mix several snp sets  
  dummy = Genome(name = 'GRCh37.75', SNPs = ['dummySRY', 'anotherSet'], SNPFilter = myFilter)


Filtering SNPs:
---------------
pyGeno allows you to select the Polymorphisms that end up into the final sequences. It supports SNPs, Inserts and Deletions.

.. code:: python

	from pyGeno.SNPFiltering import SNPFilter
	from pyGeno.SNPFiltering import SequenceSNP

	class QMax_gt_filter(SNPFilter) :
		
		def __init__(self, threshold) :
			self.threshold = threshold
			
		def filter(self, chromosome, dummySRY = None) :
			if dummySRY.Qmax_gt > self.threshold :
				#other possibilities of return are SequenceInsert(<bases>), SequenceDelete(<length>)
				return SequenceSNP(dummySRY.alt)
			return None #None means keep the reference allele
	
	persGenome = Genome(name = 'GRCh37.75_Y-Only', SNPs = 'dummySRY', SNPFilter = QMax_gt_filter(10))
	
Progress Bar:
-------------
.. code:: python

  from pyGeno.tools.ProgressBar import ProgressBar
  pg = ProgressBar(nbEpochs = 155)
  for i in range(155) :
  	pg.update(label = '%d' %i) # or simply p.update() 
  pg.close()

