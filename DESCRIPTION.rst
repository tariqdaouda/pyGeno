pyGeno: A python package for Personalized Proteogenomics
=========================================================

Short description:
------------------

With pyGeno you can do that:

.. code:: python

 from pyGeno.Genome import *
 
 #load a genome 
 ref = Genome(name = 'GRCh37.75')
 #load a gene
 gene = ref.get(Gene, name = 'TPST2')[0]
 #print the sequences of all the isoforms
 for prot in gene.get(Protein) :
  print prot.sequence

You can also do it for the **specific genomes** of your subjects:

.. code:: python

 pers = Genome(name = 'GRCh37.75', SNPs = ["RNA_S1"], SNPFilter = myFilter())

And much more: https://github.com/tariqdaouda/pyGeno

Verbose Description
--------------------

pyGeno is a python package that integrates reference sequences and annotations from Ensembl,
genomic polymorphisms from the dbSNP database and data from next-gen sequencing into an easy
to use, memory efficient and fast framework. Allowing the user to easily explore both reference 
and specific genomes and proteomes. Compared to a standalone program, pyGeno gives the user access
to the complete expressivity of python, a general programming language. It’s range of application 
therefore encompasses both short scripts and large scale genome-wide studies.

Full Documentation
------------------

The full documentation is available here_.

.. _here: http://bioinfo.iric.ca/~daoudat/pyGeno/
