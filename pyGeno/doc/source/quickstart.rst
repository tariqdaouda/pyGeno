Quickstart
==========

Quick importation
-----------------
In order to use pyGeno you first need to populate its database by importing datawraps. An importation is a one time step
and once the data has been imported, the datawrap can be discarded with no concequences to the database. 

If your goal is simply to play with pyGeno to have a feel of its magic. PyGeno is shipped with a few datawrap that require very little memory
to import.

One for importing the reference Y chromosome. This datawrap does not require any download::

	import pyGeno.bootstrap as B
	B.importHumanReference_YOnly()

There's also one that imports both the first and the Y chormosome. This tome the sequences will be imported::

	import pyGeno.bootstrap as B
	B.importHumanReference_1YOnly()

And a dummy SNP datawrap that mimics a casava's snps.txt file with only one SNP athe begining of the gene SRY::
	
	import pyGeno.bootstrap as B
	B.importDummySRY()

For more serious work pyGeno is also shipped with a datawrap for the complete human reference genome.
PyGeno will automatically download the needed data from Ensembl (~205MB). The whole importation process
may exceed an hour and requires something like 2GB of memory (according to tests)::
	
	import pyGeno.bootstrap as B
	B.importHumanReference()

That's it, you can now print the sequences of all the proteins that a gene can produce::

	from pyGeno.Genome import Genome
	from pyGeno.Gene import Gene
	from pyGeno.Protein import Protein

	#the name of the genome is defined inside the package's manifest.ini file
	ref = Genome(name = 'GRCh37.75')
	#get returns a list of elements
	gene = ref.get(Gene, name = 'SRY')[0]
	for prot in gene.get(Protein) :
		  print prot.sequence

You can see pyGeno achitecture as a graph where everything is connected to everything. For instance you can do things such as::

	gene = aProt.gene
	trans = aProt.transcript
	prot = anExon.protein
	genome = anExon.genome

Queries
-------

PyGeno allows for several kinds of queries::

	#in this case both queries will yield the same result
	myGene.get(Protein, id = "ENSID...")
	myGenome.get(Protein, id = "ENSID...")
	
	#even complex stuff
	exons = myChromosome.get(Exons, {'start >=' : x1, 'stop <' : x2})
	hlaGenes = myGenome.get(Gene, {'name like' : 'HLA'})

	sry = myGenome.get(Transcript, { "gene.name" : 'SRY' })

To know the available fields for queries, there's a "help()" function::

	Gene.help()


Faster queries
---------------

To speed up loops use iterGet()::
	
	for prot in gene.iterGet(Protein) :
	  print prot.sequence

For more speed create indexes on the fields you need the most::
	
	Gene.ensureGlobalIndex('name')


Personalized Genomes
--------------------

Personalized Genomes are a powerful feature that allow to work on the specific genomes and proteomes of your patients. You can even mix several SNPs together::
	
	from pyGeno.Genome import Genome
	#the name of the snp set is defined inside the package's manifest.ini file
	dummy = Genome(name = 'GRCh37.75', SNPs = 'dummySRY')
	#you can also define a filter (ex: a quality filter) for the SNPs
	dummy = Genome(name = 'GRCh37.75', SNPs = 'dummySRY', SNPFilter = myFilter)
	#and even mix several snp sets
	dummy = Genome(name = 'GRCh37.75', SNPs = ['dummySRY', 'anotherSet'], SNPFilter = myFilter)

pyGeno allows you to customize the Polymorphisms that end up into the final sequences. It supports SNPs, Inserts and Deletions::
	
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

