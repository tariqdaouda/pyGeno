import pyGeno.importation.Genomes as PG
import pyGeno.importation.SNPs as PS
import os

this_dir, this_filename = os.path.split(__file__)


def listDatawraps() :
	"""Lists all the datawraps pyGeno comes with"""
	l = {"Genomes" : [], "SNPs" : []}
	for f in os.listdir(os.path.join(this_dir, "bootstrap_data/genomes")) :
		if f.find(".tar.gz") > -1 :
			l["Genomes"].append(f)
	
	for f in os.listdir(os.path.join(this_dir, "bootstrap_data/SNPs")) :
		if f.find(".tar.gz") > -1 :
			l["SNPs"].append(f)

	return l

def printDatawraps() :
	"""print all available datawraps for boostraping"""
	l = listDatawraps()
	print "Available datawraps for boostraping\n"
	for k, v in l.iteritems() :
		print k
		print "~"*len(k) + "|"
		for vv in v :
			print " "*len(k) + "|" + "~~~:> " + vv
		print

def importGenome(name, batchSize = 100) :
	"""Import a genome shipped with pyGeno. Most of the datawraps only contain URLs towards data provided by third parties."""
	path = os.path.join(this_dir, "bootstrap_data", "genomes/" + name)
	PG.importGenome(path, batchSize = 100)

def importSNPs(name) :
	"""Import a SNP set shipped with pyGeno. Most of the datawraps only contain URLs towards data provided by third parties."""
	path = os.path.join(this_dir, "bootstrap_data", "genomes/" + name)
	PS.importSNPs(path, batchSize = 100)