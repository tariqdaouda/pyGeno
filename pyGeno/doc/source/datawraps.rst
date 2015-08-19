Datawraps
=========

Datawraps are used by pyGeno to import data into it's database. All reference genomes are downloaded from Ensembl, dbSNP data from dbSNP.
The :doc:`/bootstraping` module has functions to import datawraps shipped with pyGeno and also to import datawraps made available on remote location.

Importation
-----------

Here's how you import a reference genome datawrap::

	from pyGeno.importation.Genomes import *
	importGenome("my_datawrap.tar.gz")


And a SNP set datawrap::
	
	from pyGeno.importation.SNPs import *
	importSNPs("my_datawrap.tar.gz")


Creating you own datawraps
--------------------------

For polymorphims, create a file called **manifest.ini** with the following format::

	[package_infos]
	description = SNPs for testing purposes
	maintainer = Tariq Daouda
	maintainer_contact = tariq.daouda [at] umontreal
	version = 1

	[set_infos]
	species = human
	name = mySNPSET
	type = Agnostic # or CasavaSNP or dbSNPSNP
	source = Where do these snps come from?

	[snps]
	filename = snps.txt # or http://www.example.com/snps.txt or ftp://www.example.com/snps.txt if you chose not to include the file in the archive

And compress the **manifest.ini** file along sith the snps.txt (if you chose to include it and not to specify an url) into a tar.gz archive


Natively pyGeno supports dbSNP and casava(snp.txt), but it also has its own polymorphism file format (AgnosticSNP) wich is simply a tab delemited file in the following format::

	chromosomeNumber uniqueId   start        end      ref    alleles   quality       caller
	        Y          1       2655643      2655644	   T       AG        30          TopHat
	        Y          2       2655645      2655647    -       AG        28          TopHat
	        Y          3       2655648      2655650    TT      -         10          TopHat

Even tough all field are mandatory, the only ones that are critical for pyGeno to be able insert polymorphisms at the right places are: *chromosomeNumber* and *start*. All the other fields are non critical and can follow any convention you wish to apply to them, including the *end* field. You can choose the convention that suits best the query that you are planning to make on SNPs through .get(), or the way you intend to filter them using filtering objtecs.

For genomes, the manifet.ini file looks like this::

	[package_infos]
	description = Test package. This package installs only chromosome Y of mus musculus
	maintainer = Tariq Daouda
	maintainer_contact = tariq.daouda [at] umontreal
	version = GRCm38.73

	[genome]
	species = Mus_musculus
	name = GRCm38_test
	source = http://useast.ensembl.org/info/data/ftp/index.html

	[chromosome_files]
	Y = Mus_musculus.GRCm38.73.dna.chromosome.Y.fa.gz # or an url such as ftp://... or http://

	[gene_set]
	gtf = Mus_musculus.GRCm38.73_Y-only.gtf.gz # or an url such as ftp://... or http://

File URLs for refercence genomes can be found on `Ensembl: http://useast.ensembl.org/info/data/ftp/index.html`_

To learn more about datawraps and how to make your own you can have a look at :doc:`/importation`, and the Wiki_.

.. _Wiki: https://github.com/tariqdaouda/pyGeno/wiki/How-to-create-a-pyGeno-datawrap-to-import-your-data
.. _`Ensembl: http://useast.ensembl.org/info/data/ftp/index.html`: http://useast.ensembl.org/info/data/ftp/index.html
