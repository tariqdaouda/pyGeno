from pyGeno.importation.Genomes import importGenome
from pyGeno.importation.SNPs import importSNPs
import os

this_dir, this_filename = os.path.split(__file__)


def listDataWraps() :
	"""Lists all the datawraps pyGeno comes with"""
	l = []
	for f in os.listdir(os.path.join(this_dir, "bootstrap_data")) :
		if f.find(".tar.gz") > -1 :
			l.append(f)
	return l

def importHumanReference_YOnly(batchSize = 100) :
	"""Import only the Y chromosome of the Human Reference Genome (GRCh37.75). Useful for playing a bit with pyGeno.
	batchSize is the number of genes saved with each batch. Higher values mean less time wasted in io operations, but more ram needed"""
	path = os.path.join(this_dir, "bootstrap_data", "Homo_sapiens.GRCh37.75_Y-Only.tar.gz")
	importGenome(path, batchSize = 100)
 
def importHumanReference_1YOnly(batchSize = 100) :
	"""Import only the Y and first chromosomes of the Human Reference Genome (GRCh37.75). Useful for playing a bit with pyGeno.
	batchSize is the number of genes saved with each batch. Higher values mean less time wasted in io operations, but more ram needed"""
	path = os.path.join(this_dir, "bootstrap_data", "Homo_sapiens.GRCh37.75_1-Y-Only.tar.gz")
	importGenome(path, batchSize = 100)

def importDummySRY() :
	"A dummy set of SNPs for the Gene SRY on the Y chromosome."
	path = os.path.join(this_dir, "bootstrap_data", "Homo_sapiens.dummySRY.tar.gz")
	importSNPs(path)

def importHumanReference(batchSize = 100) :
	""""Import the Human Reference Genome (GRCh38.78). This may take a while, depending on the computer and 
	indexes in the database. But it's done only once.
	batchSize is the number of genes saved with each batch. Higher values mean less time wasted in io operations, but more ram needed"""
	path = os.path.join(this_dir, "bootstrap_data", "Homo_sapiens.GRCh38.78.tar.gz")
	importGenome(path, batchSize = 100)

def importMouseReference(batchSize = 100) :
	""""Import the Human Reference Genome (GRCm38.78). This may take a while, depending on the computer and 
	indexes in the database. But it's done only once.
	batchSize is the number of genes saved with each batch. Higher values mean less time wasted in io operations, but more ram needed"""
	path = os.path.join(this_dir, "bootstrap_data", "Mus_musculus.GRCm38.78.tar.gz")
	importGenome(path, batchSize = 100)

def importHumanCommonSNPs() :
	"""Import a polymorphism datawrap pyGeno comes with"""
	path = os.path.join(this_dir, "bootstrap_data", "dbSNP142_human_common_all.tar.gz")
	importSNPs(path)

def importHumanGRCh37CommonSNPs() :
	"""Import a polymorphism datawrap pyGeno comes with"""
	path = os.path.join(this_dir, "bootstrap_data", "dbSNP142_human_GRCh37_common_all.tar.gz")
	importSNPs(path)

def importGenomeDatawrap(name, batchSize = 100) :
	"""Import a genome datawrap pyGeno comes with"""
	path = os.path.join(this_dir, "bootstrap_data", name)
	importGenome(path, batchSize = 100)

def importSNPDatawrap(name, batchSize = 100) :
	"""Import a polymorphism datawrap pyGeno comes with"""
	path = os.path.join(this_dir, "bootstrap_data", name)
	importSNPs(path, batchSize = 100)
