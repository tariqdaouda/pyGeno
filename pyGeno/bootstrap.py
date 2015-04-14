import pyGeno.importation.Genomes as PG
import pyGeno.importation.SNPs as PS
from pyGeno.tools.io import printf
import os, tempfile, urllib, urllib2, json

this_dir, this_filename = os.path.split(__file__)


def listRemoteDatawraps(location = 'http://pygeno.iric.ca/_downloads/datawraps.json') :
	"""Lists all the datawraps availabe from a remote location default is 'http://pygeno.iric.ca/_downloads/datawraps.json'"""
	response = urllib2.urlopen(location)
	js = json.loads(response.read())

	return js

def printRemoteDatawraps(location) :
	"""print all available datawraps from a remote location default is 'http://pygeno.iric.ca/_downloads/datawraps.json'"""
	l = listDatawraps_url()
	printf("Available datawraps for bootstraping\n")
	for typ, dw in l.iteritems() :
		printf(typ)
		printf("~"*len(typ) + "|")
		for name in dw :
			printf(" "*len(typ) + "|" + "~~~:> " + name)
		printf('\n')

def _DW(name, url) :
	packageDir = tempfile.mkdtemp(prefix = "pyGeno_remote_")
	
	printf("~~~:>\n\tDownloading datawrap: %s..." % name)
	finalFile = os.path.normpath('%s/%s' %(packageDir, name))
	urllib.urlretrieve (url, finalFile)
	printf('\tdone.\n~~~:>')
	return finalFile

def importRemoteGenome(name, batchSize = 100) :
	"""Import a genome available from http://pygeno.iric.ca."""
	
	try :
		dw = listRemoteDatawraps()["genomes"][name]
	except AttributeError :
		raise AttributeError("There's no remote genome datawrap by the name of: '%s'" % name)

	finalFile = _DW(name, dw["url"])
	PG.importGenome(finalFile, batchSize)

def importRemoteSNPs(name) :
	"""Import a SNP set available from http://pygeno.iric.ca."""
	
	try :
		dw = listRemoteDatawraps()["snps"][name]
	except AttributeError :
		raise AttributeError("There's no remote genome datawrap by the name of: '%s'" % name)

	finalFile = _DW(name, dw["url"])
	PS.importSNPs(finalFile)

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
	"""print all available datawraps for bootstraping"""
	l = listDatawraps()
	printf("Available datawraps for boostraping\n")
	for k, v in l.iteritems() :
		printf(k)
		printf("~"*len(k) + "|")
		for vv in v :
			printf(" "*len(k) + "|" + "~~~:> " + vv)
		printf('\n')

def importGenome(name, batchSize = 100) :
	"""Import a genome shipped with pyGeno. Most of the datawraps only contain URLs towards data provided by third parties."""
	path = os.path.join(this_dir, "bootstrap_data", "genomes/" + name)
	PG.importGenome(path, batchSize)

def importSNPs(name) :
	"""Import a SNP set shipped with pyGeno. Most of the datawraps only contain URLs towards data provided by third parties."""
	path = os.path.join(this_dir, "bootstrap_data", "SNPs/" + name)
	PS.importSNPs(path)