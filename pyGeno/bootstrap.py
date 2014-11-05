from .importation.Genomes import importGenome
from .importation.SNPs import importSNPs
import os
this_dir, this_filename = os.path.split(__file__)

def importHumanReference() :
	path = os.path.join(this_dir, "bootstrap_data", "GRCh37.75.tar.gz")
	importGenome(path)

def importHumanReference_YOnly() :
	path = os.path.join(this_dir, "bootstrap_data", "GRCh37.75_Y-Only.tar.gz")
	importGenome(path)

def importDummySRY() :
	path = os.path.join(this_dir, "bootstrap_data", "dummySRY.tar.gz")
	importSNPs(path)
