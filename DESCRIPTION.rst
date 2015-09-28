pyGeno: A python package for Precision Medicine, Personalized Genomics and Proteomics
=====================================================================================

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

pyGeno is a personal bioinformatic database that runs directly into python, on your laptop and does not depend
upon any REST API. pyGeno is here to make extracting data such as gene sequences a breeze, and is designed to
be able cope with huge queries. The most exciting feature of pyGeno, is that it allows to work with seamlessly with both reference and **Presonalized Genomes**.

Personalized Genomes, are custom genomes that you create by combining a reference genome, sets of polymorphims and an optional filter.
pyGeno will take care of applying the filter and inserting the polymorphisms at their right place, so you get
direct access to the DNA and Protein sequences of your patients.

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
