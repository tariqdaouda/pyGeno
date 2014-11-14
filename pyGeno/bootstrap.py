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
	"""Importing only the Y chromosome of the Human Reference Genome. Useful for playing a bit with pyGeno.
	batchSize is the number of genes saved with each batch. Higher values mean less time wasted in io operations, but more ram needed"""
	path = os.path.join(this_dir, "bootstrap_data", "Homo_sapiens.GRCh37.75_Y-Only.tar.gz")
	importGenome(path, batchSize = 100)
 
def importHumanReference_1YOnly(batchSize = 100) :
	"""Importing only the Y and first chromosomes of the Human Reference Genome. Useful for playing a bit with pyGeno.
	batchSize is the number of genes saved with each batch. Higher values mean less time wasted in io operations, but more ram needed"""
	path = os.path.join(this_dir, "bootstrap_data", "Homo_sapiens.GRCh37.75_1-Y-Only.tar.gz")
	importGenome(path, batchSize = 100)

def importDummySRY() :
	"A dummy set of SNPs for the Gene SRY on the Y chromosome."
	path = os.path.join(this_dir, "bootstrap_data", "Homo_sapiens.dummySRY.tar.gz")
	importSNPs(path)

def importHumanReference(batchSize = 100) :
	""""Importes the Human Reference Genome. This may take a while, depending on the computers and 
	indexes in the database. But it's done only once.
	batchSize is the number of genes saved with each batch. Higher values mean less time wasted in io operations, but more ram needed"""
	path = os.path.join(this_dir, "bootstrap_data", "Homo_sapiens.GRCh37.75.tar.gz")
	importGenome(path, batchSize = 100)
