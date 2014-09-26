pyGeno: A python package for Personalized Proteogenomics
========================================================

Bootstraping:
-------------
There are two packages shipped with pyGeno, the reference genome GRCh37.75 and a dummy packages that mimics a casava's snps.txt with one SNP at begining of the gene SRY. 
You can bootstrap pyGeno with the following

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

For more details about how packages are built you can have a look inside the folder bootstrap_data. There are two packages for you there.

Instanciating a genome:
-----------------------
.. code:: python
	
	from pyGeno.Genome import Genome
	ref = Genome(name = 'GRCh37.75')

Printing all the proteins of a gene:
-----------------------------------
.. code:: python

  from pyGeno.Genome import Genome
  from pyGeno.Gene import Gene
  from pyGeno.Protein import Protein
  ref = Genome(name = 'GRCh37.75')
  gene = ref.get(Gene, name = 'TPST2')[0]
  for prot in gene.get(Protein) :
  	print prot.sequence

Making queries:
--------------
pyGeno's get function uses the expressivity of rabaDB

These are all possible query formats

.. code:: python
  ref.get(Gene, name = "SRY")
  ref.get(Gene, { "name like" : "HLA"})
  chormosome.get(Exon, { "start >" : 12000, "end <" : 12300 })

Making queries get vs iterGet:
-----------------------------
iterGet is a faster version of get that returns an iterator instead of a list.


Creating indexes to speed up queries:
------------------------------------
.. code:: python

  from pyGeno.Gene import Gene
  #creating an index on gene names is necessary
  Gene.ensureGobalIndex('name')
  
Loading a genome with SNPs:
---------------------------
.. code:: python
  
  from pyGeno.Genome import Genome
  #the name of the snp set is defined inside the package
  dummy = Genome(name = 'GRCh37.75', SNPs = 'dummySRY')
  #you can also define a filter
  dummy = Genome(name = 'GRCh37.75', SNPs = 'dummySRY', SNPFilter = myFilter)
  #and mix several snp sets  
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

