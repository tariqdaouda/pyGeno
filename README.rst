pyGeno: Personalized Genomics and Protemomics
=============================================

.. image:: http://bioinfo.iric.ca/~daoudat/pyGeno/_static/logo.png
   :alt: pyGeno's logo

pyGeno is developed by `Tariq Daouda`_ at the *Institute for Research in Immunology and Cancer* (IRIC_), its logo is the work of the freelance designer `Sawssan Kaddoura`_.
For the latest news about pyGeno, you can follow me on twitter `@tariqdaouda`_.

.. _Tariq Daouda: http://wwww.tariqdaouda.com
.. _IRIC: http://www.iric.ca
.. _Sawssan Kaddoura: http://www.sawssankaddoura.com


Click here for The `full documentation`_.

.. _full documentation: http://pygeno.iric.ca/

For the latest news about pyGeno, you can follow me on twitter `@tariqdaouda`_.

.. _@tariqdaouda: https://www.twitter.com/tariqdaouda

Installation:
-------------

.. code:: python
	
	pip install pyGeno #for the latest stable version

For the latest developments:

.. code:: shell

	git clone https://github.com/tariqdaouda/pyGeno.git
	cd pyGeno
	python setup.py develop

A brief introduction
--------------------

pyGeno is a personal bioinformatic database that runs directly into python, on your laptop and does not depend
upon any REST API. pyGeno is here to make extracting data such as gene sequences a breeze, and is designed to
be able cope with huge queries.

.. code:: python

	from pyGeno.Genome import *
	
	g = Genome(name = "GRCh37.75")
	prot = g.get(Protein, id = 'ENSP00000438917')[0]
	#print the protein sequence
	print prot.sequence
	#print the protein's gene biotype
	print prot.gene.biotype
	#print protein's transcript sequence
	print protein.transcript.sequence
	
	...
	
	#fancy queries
	for exons in g.get(Exons, {"CDS_start >": x1, "CDS_end <=" : x2, "chromosome.number" : "22"}) :
		#print the exon's coding sequence
		print exon.CDS
		#print the exon's transcript sequence
		print exon.transcript.sequence
	
	...
	
	#You can do the same for your subject specific genomes
	#by combining a reference genome with polymorphisms 
	g = Genome(name = "GRCh37.75", SNPs = ["STY21_RNA"], SNPFilter = MyFilter())

And if you ever get lost, there's an online **help()** function for each object type:

.. code:: python

	from pyGeno.Genome import *
	
	print Exon.help()

Should output:

.. code::
	
	Available fields for Exon: CDS_start, end, chromosome, CDS_length, frame, number, CDS_end, start, genome, length, protein, gene, transcript, id, strand

Batteries included (bootstraping):
---------------------------------

pyGeno's database is populated by importing data wraps.
pyGeno comes with a few data wraps, to get the list you can use:

.. code:: python
	
	import pyGeno.bootstrap as B
	B.printDatawraps()

.. code::

	Available datawraps for boostraping
	
	SNPs
	~~~~|
	    |~~~:> Homo_sapiens_agnostic.dummySRY.tar.gz
	    |~~~:> Homo_sapiens.dummySRY_casava.tar.gz
	    |~~~:> dbSNP142_human_GRCh37_common_all.tar.gz
	    |~~~:> dbSNP142_human_common_all.tar.gz
	
	
	Genomes
	~~~~~~~|
	       |~~~:> Homo_sapiens.GRCh37.75.tar.gz
	       |~~~:> Homo_sapiens.GRCh37.75_Y-Only.tar.gz
	       |~~~:> Homo_sapiens.GRCh38.78.tar.gz
	       |~~~:> Mus_musculus.GRCm38.78.tar.gz

To get a list of remote datawraps that pyGeno can download for you, do:

.. code:: python

	B.printRemoteDatawraps()

Importing whole genomes is a demanding process that take more than an hour and requires (according to tests) 
at least 3GB of memory. Depending on your configuration, more might be required.

That being said importating a data wrap is a one time operation and once the importation is complete the datawrap
can be discarded without consequences.

The bootstrap module also has some handy functions for importing built-in packages.

Some of them just for playing around with pyGeno (**Fast importation** and **Small memory requirements**):

.. code:: python
	
	import pyGeno.bootstrap as B

	#Imports only the Y chromosome from the human reference genome GRCh37.75
	#Very fast, requires even less memory. No download required.
	B.importGenome("Homo_sapiens.GRCh37.75_Y-Only.tar.gz")
	
	#A dummy datawrap for humans SNPs and Indels in pyGeno's AgnosticSNP  format. 
	# This one has one SNP at the begining of the gene SRY
	B.importSNPs("Homo_sapiens.dummySRY_casava.tar.gz")

And for more **Serious Work**, the whole reference genome.

.. code:: python

	#Downloads the whole genome (205MB, sequences + annotations), may take an hour or more.
	B.importGenome("Homo_sapiens.GRCh38.78.tar.gz")
	
Importing a custom datawrap:
--------------------------

.. code:: python

  from pyGeno.importation.Genomes import *
  importGenome('GRCh37.75.tar.gz')

To import a patient's specific polymorphisms

.. code:: python

  from pyGeno.importation.SNPs import *
  importSNPs('patient1.tar.gz')

For a list of available datawraps available for download, please have a look here_.

You can easily make your own datawraps with any tar.gz compressor.
For more details on how datawraps are made you can check wiki_ or have a look inside the folder bootstrap_data.

.. _here: http://pygeno.iric.ca/datawraps.html
.. _wiki: https://github.com/tariqdaouda/pyGeno/wiki/How-to-create-a-pyGeno-friendly-package-to-import-your-data%3F

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

Or simply:

.. code:: python

  from pyGeno.Genome import *

then:

.. code:: python

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

Find in sequences:
------------------

Internally pyGeno uses a binary representation for nucleotides and amino acids to deal with polymorphisms. 
For example,both "AGC" and "ATG" will match the following sequence "...AT/GCCG...".

.. code:: python

	#returns the position of the first occurence
	transcript.find("AT/GCCG")
	#returns the positions of all occurences
	transcript.findAll("AT/GCCG")
	
	#similarly, you can also do
	transcript.findIncDNA("AT/GCCG")
	transcript.findAllIncDNA("AT/GCCG")
	transcript.findInUTR3("AT/GCCG")
	transcript.findAllInUTR3("AT/GCCG")
	transcript.findInUTR5("AT/GCCG")
	transcript.findAllInUTR5("AT/GCCG")
	
	#same for proteins
	protein.find("DEV/RDEM")
	protein.findAll("DEV/RDEM")
	
	#and for exons
	exon.find("AT/GCCG")
	exon.findAll("AT/GCCG")
	exon.findInCDS("AT/GCCG")
	exon.findAllInCDS("AT/GCCG")
	...
	
Creating a Personalized Genome:
-------------------------------
Personalized Genomes are a powerful feature that allow to work on the specific genomes and proteomes of your patients.
You can even mix several SNPs together.

.. code:: python
  
  from pyGeno.Genome import Genome
  #the name of the snp set is defined inside the package's manifest.ini file
  dummy = Genome(name = 'GRCh37.75', SNPs = 'dummySRY')
  #you can also define a filter (ex: a quality filter) for the SNPs
  dummy = Genome(name = 'GRCh37.75', SNPs = 'dummySRY', SNPFilter = myFilter())
  #and even mix several snp sets  
  dummy = Genome(name = 'GRCh37.75', SNPs = ['dummySRY', 'anotherSet'], SNPFilter = myFilter())


Filtering SNPs:
---------------
pyGeno allows you to select the Polymorphisms that end up into the final sequences. It supports SNPs, Inserts and Deletions.

.. code:: python

	from pyGeno.SNPFiltering import SNPFilter, SequenceSNP

	class QMax_gt_filter(SNPFilter) :
		
		def __init__(self, threshold) :
			self.threshold = threshold
			
		def filter(self, chromosome, dummySRY = None) :
			if dummySRY.Qmax_gt > self.threshold :
				#other possibilities of return are SequenceInsert(<bases>), SequenceDelete(<length>)
				return SequenceSNP(dummySRY.alt)
			return None #None means keep the reference allele
	
	persGenome = Genome(name = 'GRCh37.75_Y-Only', SNPs = 'dummySRY', SNPFilter = QMax_gt_filter(10))

Getting an arbitrary sequence:
------------------------------
You can ask for any sequence of any chromosome:

.. code:: python
	
	chr12 = myGenome.get(Chromosome, number = "12")[0]
	print chr12.sequence[x1:x2]
	# for the reference sequence
  	print chr12.refSequence[x1:x2]
	
Progress Bar:
-------------
.. code:: python

  from pyGeno.tools.ProgressBar import ProgressBar
  pg = ProgressBar(nbEpochs = 155)
  for i in range(155) :
  	pg.update(label = '%d' %i) # or simply p.update() 
  pg.close()

