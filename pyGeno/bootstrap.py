import pyGeno.importation.Genomes as PG
# import pyGeno.importation.SNPs as PS
from pyGeno.tools.io import printf
import os, tempfile, json
import urllib.request, urllib.parse, urllib.error, urllib.request, urllib.error, urllib.parse
import pyGeno.configuration as conf

this_dir, this_filename = os.path.split(__file__)

def _DW(name, url) :
    packageDir = tempfile.mkdtemp(prefix = "pyGeno_remote_")
    printf("~~~:>\n\tDownloading datawrap: %s..." % name)
    finalFile = os.path.normpath('%s/%s' %(packageDir, name))
    urllib.request.urlretrieve (url, finalFile)
    printf('\tdone.\n~~~:>')
    return finalFile

def list_datawraps() :
    """Lists all the datawraps pyGeno comes with"""
    l = {"Genomes" : [], "SNPs" : []}
    for f in os.listdir(os.path.join(this_dir, "bootstrap_data/genomes")) :
        if f.find(".tar.gz") > -1 :
            l["Genomes"].append(f)
    
    for f in os.listdir(os.path.join(this_dir, "bootstrap_data/SNPs")) :
        if f.find(".tar.gz") > -1 :
            l["SNPs"].append(f)

    return l

def print_datawraps() :
    """print all available datawraps for bootstraping"""
    l = listDatawraps()
    printf("Available datawraps for boostraping\n")
    for k, v in l.items() :
        printf(k)
        printf("~"*len(k) + "|")
        for vv in v :
            printf(" "*len(k) + "|" + "~~~:> " + vv)
        printf('\n')

def import_genome(name, batchSize = 100) :
    """Import a genome shipped with pyGeno. Most of the datawraps only contain URLs towards data provided by third parties."""
    path = os.path.join(this_dir, "bootstrap_data", "genomes/" + name)
    PG.import_genome(path, batchSize)

def import_SNPs(name) :
    """Import a SNP set shipped with pyGeno. Most of the datawraps only contain URLs towards data provided by third parties."""
    path = os.path.join(this_dir, "bootstrap_data", "SNPs/" + name)
    PS.import_SNPs(path)
