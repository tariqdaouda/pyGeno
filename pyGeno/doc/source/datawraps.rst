Datawraps
=========

Datawraps are used by pyGeno to import data into it's database. All reference genomes are downloaded from Ensembl, dbSNP data from dbSNP.

Here's how you import a reference genome datawrap::

	from pyGeno.importation.Genomes import *
	importGenome("my_datawrap.tar.gz")


And a SNP set datawrap::
	
	from pyGeno.importation.SNPs import *
	importSNPs("my_datawrap.tar.gz")

Natively pyGeno supports dbSNP and casava(snp.txt), but it also has its own polymorphism file format wich is simply a tab delemited file in the following format::

	chromosomeNumber	start	end	ref	alleles	quality	caller
	Y	2655643	2655644	T	AG	30	TopHat
	Y	2655645	2655647	-	AG	30	TopHat
	Y	2655648	2655650	TT	-	30	TopHat

To learn more about datawraps and how to make your own you can have a look at :doc:`/importation`, and the Wiki_.

.. _Wiki: https://github.com/tariqdaouda/pyGeno/wiki/How-to-create-a-pyGeno-datawrap-to-import-your-data

Human reference genomes
------------------------

* :download:`GRCh37.75 <./datawraps/Homo_sapiens.GRCh37.75.tar.gz>`

* :download:`GRCh38.78 <./datawraps/Homo_sapiens.GRCh38.78.tar.gz>`

Mouse reference genomes
------------------------

* :download:`GRCm38.78 <./datawraps/Mus_musculus.GRCm38.78.tar.gz>`

dbSNP
-------

**dbSNP regularly drops the support of previous versions. You might need to update the URLs in the manifest.ini file of the datawrap. If you do so, please send it to me an email so I can update this page**

* :download:`Human common SNPs for build 142 <./datawraps/dbSNP142_human_common_all.tar.gz>`
