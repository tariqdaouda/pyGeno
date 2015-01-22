Importation
===========
pyGeno's database is populated by importing tar.gz compressed archives called *datawraps*. An importation is a one time step and once a datawrap has been imported it can be discarded with no concequences to the database. 

Here's how you import a reference genome datawrap::
	
	from pyGeno.importation.Genomes import *
	importGenome("my_genome_datawrap.tar.gz")


And a SNP set datawrap::
	
	from pyGeno.importation.SNPs import *
	importSNPs("my_snp_datawrap.tar.gz")

pyGeno comes with a few datawraps that you can quickly import using the :doc:`/bootstraping` module.

You can find a list of datawraps to import here: :doc:`/datawraps`

You can also easily create your own by simply putting a bunch of URLs into a *manifest.ini* file and compressing int into a *tar.gz archive* (as explained below or on the Wiki_).

.. _Wiki: https://github.com/tariqdaouda/pyGeno/wiki/How-to-create-a-pyGeno-datawrap-to-import-your-data

Genomes
-------
.. automodule:: importation.Genomes
   :members:

Polymorphisms (SNPs and Indels)
-------------------------------

.. automodule:: importation.SNPs
   :members:
