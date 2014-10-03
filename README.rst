pyGeno: A python package for Personalized Proteogenomics
========================================================
.. code:: python
	
	pip install pyGeno #for the latest stable version


Bootstraping:
-------------
There are two data wraps shipped with pyGeno, the reference genome GRCh37.75 and a dummy data wrap that mimics a casava's snps.txt with one SNP at the begining of the gene SRY. 
You can bootstrap pyGeno with the following:

.. code:: python
	
	import pyGeno.bootstrap as B
	
	B.importHumanReference()
	B.importDummySRY()

Importing a data wrap:
----------------------

.. code:: python

  from pyGeno.importation.Genomes import *
  importGenome('GRCh37.75.tar.gz')

To import a patient's specific polymorphisms

.. code:: python

  from pyGeno.importation.SNPs import *
  importSNPs('patient1.tar.gz')

You can easily make your own data wraps with any tar.gz compressor.
For more details on how data wraps are made you can have a look inside the folder bootstrap_data. There are two there for you.

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
  Gene.ensureGobalIndex('name')
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
For an example of how to define your own filters you can have a look at the function defaultSNPFilter in SNP.py

Progress Bar:
-------------
.. code:: python

  from pyGeno.tools.ProgressBar import ProgressBar
  pg = ProgressBar(nbEpochs = 155)
  for i in range(155) :
  	p.update(label = '%d' %i) # or simply p.update() 
  p.close()

