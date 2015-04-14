import urllib, shutil

from ConfigParser import SafeConfigParser
import pyGeno.configuration as conf
from pyGeno.SNP import *
from pyGeno.tools.ProgressBar import ProgressBar
from pyGeno.tools.io import printf
from Genomes import _decompressPackage, _getFile

from pyGeno.tools.parsers.CasavaTools import SNPsTxtFile
from pyGeno.tools.parsers.VCFTools import VCFFile
from pyGeno.tools.parsers.CSVTools import CSVFile

def importSNPs(packageFile) :
	"""The big wrapper, this function should detect the SNP type by the package manifest and then launch the corresponding function.
	Here's an example of a SNP manifest file for Casava SNPs::

		[package_infos]
		description = Casava SNPs for testing purposes
		maintainer = Tariq Daouda
		maintainer_contact = tariq.daouda [at] umontreal
		version = 1

		[set_infos]
		species = human
		name = dummySRY
		type = Agnostic
		source = my place at IRIC

		[snps]
		filename = snps.txt # as with genomes you can either include de file at the root of the package or specify an URL from where it must be downloaded
	"""
	printf("Importing polymorphism set: %s... (This may take a while)" % packageFile)
	
	packageDir = _decompressPackage(packageFile)

	parser = SafeConfigParser()
	parser.read(os.path.normpath(packageDir+'/manifest.ini'))
	packageInfos = parser.items('package_infos')

	setName = parser.get('set_infos', 'name')
	typ = parser.get('set_infos', 'type')+'SNP'
	species = parser.get('set_infos', 'species').lower()
	genomeSource = parser.get('set_infos', 'source')
	snpsFileTmp = parser.get('snps', 'filename').strip()
	snpsFile = _getFile(parser.get('snps', 'filename'), packageDir)
	
	try :
		SMaster = SNPMaster(setName = setName)
	except KeyError :
		if typ.lower() == 'casavasnp' :
			return _importSNPs_CasavaSNP(setName, species, genomeSource, snpsFile)
		elif typ.lower() == 'dbsnpsnp' :
			return _importSNPs_dbSNPSNP(setName, species, genomeSource, snpsFile)
		elif typ.lower() == 'topHatsnp' :
			return _importSNPs_TopHatSNP(setName, species, genomeSource, snpsFile)
		elif typ.lower() == 'agnosticsnp' :
			return _importSNPs_AgnosticSNP(setName, species, genomeSource, snpsFile)
		else :
			raise FutureWarning('Unknown SNP type in manifest %s' % typ)
	else :
		raise KeyError("There's already a SNP set by the name %s. Use deleteSNPs() to remove it first" %setName)
	
	shutil.rmtree(packageDir)

def deleteSNPs(setName) :
	"""deletes a set of polymorphisms"""
	con = conf.db
	try :
		SMaster = SNPMaster(setName = setName)
		con.beginTransaction()
		SNPType = SMaster.SNPType
		con.delete(SNPType, 'setName = ?', (setName,))
		SMaster.delete()
		con.endTransaction()
	except KeyError :
		raise KeyError("Can't delete the setName %s because i can't find it in SNPMaster, maybe there's not set by that name" % setName)
		#~ printf("can't delete the setName %s because i can't find it in SNPMaster, maybe there's no set by that name" % setName)
		return False
	return True

def _importSNPs_AgnosticSNP(setName, species, genomeSource, snpsFile) :
	"This function will also create an index on start->chromosomeNumber->setName. Warning : pyGeno wil interpret all positions as 0 based"
	printf('importing SNP set %s for species %s...' % (setName, species))

	snpData = CSVFile()
	snpData.parse(snpsFile, separator = "\t")

	AgnosticSNP.dropIndex(('start', 'chromosomeNumber', 'setName'))
	conf.db.beginTransaction()
	
	pBar = ProgressBar(len(snpData))
	pLabel = ''
	currChrNumber = None
	for snpEntry in snpData :
		tmpChr = snpEntry['chromosomeNumber']
		if tmpChr != currChrNumber :
			currChrNumber = tmpChr
			pLabel = 'Chr %s...' % currChrNumber

		snp = AgnosticSNP()
		snp.species = species
		snp.setName = setName
		for f in snp.getFields() :
			try :
				setattr(snp, f, snpEntry[f])
			except KeyError :
				if f != 'species' and f != 'setName' :
					printf("Warning filetype as no key %s", f)
		snp.start = int(snp.start)
		snp.end = int(snp.end)
		snp.save()
		pBar.update(label = pLabel)

	pBar.close()
	
	snpMaster = SNPMaster()
	snpMaster.set(setName = setName, SNPType = 'AgnosticSNP', species = species)
	snpMaster.save()

	printf('saving...')
	conf.db.endTransaction()
	printf('creating indexes...')
	CasavaSNP.ensureGlobalIndex(('start', 'chromosomeNumber', 'setName'))
	printf('importation of SNP set %s for species %s done.' %(setName, species))
	
	return True

def _importSNPs_CasavaSNP(setName, species, genomeSource, snpsFile) :
	"This function will also create an index on start->chromosomeNumber->setName. Warning : pyGeno positions are 0 based"
	printf('importing SNP set %s for species %s...' % (setName, species))

	snpData = SNPsTxtFile(snpsFile)
	
	CasavaSNP.dropIndex(('start', 'chromosomeNumber', 'setName'))
	conf.db.beginTransaction()
	
	pBar = ProgressBar(len(snpData))
	pLabel = ''
	currChrNumber = None
	for snpEntry in snpData :
		tmpChr = snpEntry['chromosomeNumber']
		if tmpChr != currChrNumber :
			currChrNumber = tmpChr
			pLabel = 'Chr %s...' % currChrNumber

		snp = CasavaSNP()
		snp.species = species
		snp.setName = setName
		
		for f in snp.getFields() :
			try :
				setattr(snp, f, snpEntry[f])
			except KeyError :
				if f != 'species' and f != 'setName' :
					printf("Warning filetype as no key %s", f)
		snp.start -= 1
		snp.end -= 1
		snp.save()
		pBar.update(label = pLabel)

	pBar.close()
	
	snpMaster = SNPMaster()
	snpMaster.set(setName = setName, SNPType = 'CasavaSNP', species = species)
	snpMaster.save()

	printf('saving...')
	conf.db.endTransaction()
	printf('creating indexes...')
	CasavaSNP.ensureGlobalIndex(('start', 'chromosomeNumber', 'setName'))
	printf('importation of SNP set %s for species %s done.' %(setName, species))
	
	return True

def _importSNPs_dbSNPSNP(setName, species, genomeSource, snpsFile) :
	"This function will also create an index on start->chromosomeNumber->setName. Warning : pyGeno positions are 0 based"
	snpData = VCFFile(snpsFile, gziped = True, stream = True)
	dbSNPSNP.dropIndex(('start', 'chromosomeNumber', 'setName'))
	conf.db.beginTransaction()
	pBar = ProgressBar()
	pLabel = ''
	for snpEntry in snpData :
		pBar.update(label = 'Chr %s, %s...' %  (snpEntry['#CHROM'], snpEntry['ID']))
		
		snp = dbSNPSNP()
		for f in snp.getFields() :
			try :
				setattr(snp, f, snpEntry[f])
			except KeyError :
				pass
		snp.chromosomeNumber = snpEntry['#CHROM']
		snp.species = species
		snp.setName = setName
		snp.start = snpEntry['POS']-1
		snp.alt = snpEntry['ALT']
		snp.ref = snpEntry['REF']
		snp.end = snp.start+len(snp.alt)
		snp.save()
	
	pBar.close()
	
	snpMaster = SNPMaster()
	snpMaster.set(setName = setName, SNPType = 'dbSNPSNP', species = species)
	snpMaster.save()
	
	printf('saving...')
	conf.db.endTransaction()
	printf('creating indexes...')
	dbSNPSNP.ensureGlobalIndex(('start', 'chromosomeNumber', 'setName'))
	printf('importation of SNP set %s for species %s done.' %(setName, species))

	return True
	
def _importSNPs_TopHatSNP(setName, species, genomeSource, snpsFile) :
	raise FutureWarning('Not implemented yet')
