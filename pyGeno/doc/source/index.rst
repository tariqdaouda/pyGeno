.. pyGeno documentation master file, created by
   sphinx-quickstart on Thu Nov  6 16:45:34 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyGeno: Ensembl, Next-Gen Sequencing, dbSNP together on your laptop
===================================================================

pyGeno's `lair is on Github`_.

.. _lair is on Github: http://www.github.com/tariqdaouda/pyGeno


A compelling example:
----------------------

.. code::

	from pyGeno.Genome import *

	g = Genome(name = "GRCh37.75")
	prot = g.get(Protein, id = 'ENSP00000438917')[0]
	print prot.sequence
	print prot.gene.biotype
	print protein.transcript.sequence

	...

	#fancy queries
	for exons in g.get(Exons, {"CDS_start >": x1, "CDS_end <=" : x2, "chromosome.number" : "22"}) :
		    print exon.CDS

	...

	#and you can do the same for your subject specific genomes
	g = Genome(name = "GRCh37.75", SNPs = ["STY21_RNA"], SNPFilter = MyFilter())


Verbose Introduction
---------------------

pyGeno integrates:

* **Reference sequences** and annotations from **Ensembl**

* Genomic polymorphisms from the **dbSNP** database

* SNPs from **next-gen** sequencing

pyGeno was designed to be:

* Fast to install. It has no dependencies but its own backend: `rabaDB`_.
* Fast to run and memory efficient, so you can use it on your laptop.
* Fast to use. No queries to foreign APIs all the data rests on your computer, so it is readily accessible when you need it.
* Fast to learn. One sigle function **get()** can do the job of several other tools at once. 

It also comes with:

* Parsers for: FASTA, FASTQ, GTF, VCF, CSV.
* Useful tools for traduction etc...
* Optimised genome indexation with *Segment Trees*.
* A funky *Progress Bar*.

But in my opinion the coolest thing about pyGeno is that it allows to create **personalized genomes**. Genomes that you design yourself by combining reference genomes and SNP sets derived next-gen sequencing or dbSNP.

Give it a try and let me know what you think!

Cheers,

pyGeno is developed by `Tariq Daouda`_ at the *Institute for Research in Immunology and Cancer* (IRIC_), its logo is the work of the freelance designer `Sawssan Kaddoura`_.
For the latest news about pyGeno, you can follow me on twitter `@tariqdaouda`_.

.. _rabaDB: https://github.com/tariqdaouda/rabaDB
.. _Tariq Daouda: http://www.tariqdaouda.com
.. _IRIC: http://www.iric.ca
.. _Sawssan Kaddoura: http://www.sawssankaddoura.com
.. _@tariqdaouda: https://www.twitter.com/tariqdaouda

Contents:
----------

.. toctree::
   :maxdepth: 2
   
   publications
   quickstart
   installation
   bootstraping
   datawraps
   importation
   objects
   snp_filter
   tools
   parsers

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

