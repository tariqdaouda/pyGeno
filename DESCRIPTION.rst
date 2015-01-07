pyGeno: A python package for Personalized Medicine and Proteogenomics
=====================================================================

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

pyGeno is mainly intended for personalized medicine, but it can serve many purposes. It integrates reference sequences and 
annotations, genomic polymorphisms from the dbSNP database 
and data from next-gen sequencing into an easy to use, 
memory- efficient and fast framework, therefore allowing 
the user to easily explore reference and custom genomes and 
proteomes. pyGeno can be used for both short scripts and large scale genome-wide studies.

Full Documentation
------------------

The full documentation is available here_

.. _here: http://pygeno.iric.ca/

For the latest news about pyGeno, you can follow me on twitter `@tariqdaouda`_.

.. _@tariqdaouda: https://www.twitter.com/tariqdaouda
