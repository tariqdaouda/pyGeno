pyGeno: A python package for Personalized Genomics and Proteomics
=================================================================

Short description:
------------------

pyGeno is developed by `Tariq Daouda`_ at the *Institute for Research in Immunology and Cancer* (IRIC_).

.. _Tariq Daouda: http://www.tariqdaouda.com
.. _IRIC: http://www.iric.ca

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

pyGeno is a python package that allows you to easily combine Reference Genomes and sets of Polymorphisms together to create personalized genomes. Personalized genomes can be used to work directly on the genomes of you subjects and be translated into Personalized Proteomes, 

Multiple sets of of polymorphisms can also be combined together to leverage their independent benefits ex: 

RNA-seq and DNA-seq for the same individual to improve the coverage
RNA-seq of an individual + dbSNP for validation
Combine the results of RNA-seq of several individual to create a genome only containing the common polymorphisms
pyGeno is also a personal database that give you access to all the information provided by Ensembl (for both Reference and Personalized Genomes) without the need of queries to distant HTTP APIs. Allowing for much faster and reliable genome wide study pipelines.

It also comes with parsers for several file types and various other useful tools.

Full Documentation
------------------

The full documentation is available here_

.. _here: http://pygeno.iric.ca/

If you like pyGeno, please let me know.
For the latest news, you can follow me on twitter `@tariqdaouda`_.

.. _@tariqdaouda: https://www.twitter.com/tariqdaouda
