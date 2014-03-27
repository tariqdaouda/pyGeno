import os, glob, gzip, tarfile, shutil, time, sys, gc, cPickle, tempfile
from ConfigParser import SafeConfigParser

import configuration as conf

from SNP import *

from pyGeno.tools.ProgressBar import ProgressBar
from pyGeno.tools.io import printf
from Genomes import _decompressPackage

from pyGeno.tools.CasavaTools import CasavaTools

def importSNPs(packageFile) :
	"""The big wrapper, this function should detect the SNP type by the package manifest and then launch the corresponding function"""
	packageDir = _decompressPackage(packageFile)

	parser = SafeConfigParser()
	parser.read(os.path.normpath(packageDir+'/manifest.ini'))
	packageInfos = parser.items('package_infos')

	setName = parser.get('set', 'name')
	typ = setName = parser.get('set', 'type')
	specie = parser.get('set', 'specie').lower()
	genomeSource = parser.get('set', 'source')
	snpsTxtFile = os.path.normpath('%s/%s' %(packageDir, parser.get('snps', 'casava_snps_txt')))
	
	if typ == 'CasavaSNP' :
		_importSNPs_CasavaSNP(setName, specie, genomeSource, snpsTxtFile)
	elif typ == 'dbSNPSNP' :
		raise FutureWarning('Not implemented yet')
	elif typ == 'TopHatSNP' :
		raise FutureWarning('Not implemented yet')
	else :
		raise FutureWarning('Unknown SNP type in manifest %s' % typ)

def deleteSNPs(setName) :
	con = conf.db
	try :
		SMaster = SNPMaster(setName = setName)
		con.beginTransaction()
		SNPType = SMaster.SNPType
		con.delete(SNPType, 'setName = ?', (setName,))
		SMaster.delete()
		con.endTransaction()
	except KeyError :
		printf("can't delete the setName %s because i can't find it in SNPMaster, maybe there's not set by that name" % setName)
		return

def _importSNPs_CasavaSNP(setName, specie, genomeSource, snpsTxtFile) :

	printf('importing SNP set %s for specie %s...' % (setName, specie))

	conf.db.beginTransaction()
	
	snpData = SNPsTxtFile(snpsTxtFile)
	
	CasavaSNP.dropIndex('setName')
	CasavaSNP.dropIndex('start')

	pBar = ProgressBar(len(lines))
	pLabel = ''
	currChrNumber = None
	for snpEntry in snpData :
		tmpChr = snpEntry['chromosomeNumber']
		if tmpChr != currChrNumber :
			currChrNumber = tmpChr
			pLabel = 'Chr %s...' % currChrNumber

		snp = CasavaSNP()
		snp.chromosomeNumber = currChrNumber
		snp.specie = specie
		snp.setName = setName
		#first column: chro, second first of range (identical to second column)
		snp.start = snpEntry['start']
		snp.end = snpEntry['end']
		snp.bcalls_used = snpEntry['bcalls_used']
		snp.bcalls_filt = snpEntry['bcalls_filt']
		snp.ref = snpEntry['ref']
		snp.QSNP = snpEntry['QSNP']
		snp.alleles = snpEntry['alleles']
		snp.Qmax_gt = snpEntry['Qmax_gt']
		snp.max_gt_poly_site = snpEntry['max_gt_poly_site']
		snp.Qmax_gt_poly_site = snpEntry['Qmax_gt_poly_site']
		snp.A_used = snpEntry['A_used']
		snp.C_used = snpEntry['C_used']
		snp.G_used = snpEntry['G_used']
		snp.T_used = snpEntry['T_used']
		snp.save()
		pBar.update(label = pLabel)

	snpMaster = SNPMaster()
	snpMaster.set(setName = setName, SNPType = 'CasavaSNP', specie = specie)
	snpMaster.save()

	printf('saving...')
	conf.db.endTransaction()
	printf('creating indexes...')
	CasavaSNP.enureIndex('setName')
	CasavaSNP.enureIndex('start')
	printf('importation of SNP set %s for specie %s done.' %(setName, specie))

def _importSNPs_dbSNPSNP(setName, specie, genomeSource, snpsTxtFile) :
	raise FutureWarning('Not implemented yet')

def _importSNPs_TopHatSNP(setName, specie, genomeSource, snpsTxtFile) :
	raise FutureWarning('Not implemented yet')
	
if __name__ == "__main__" :
	print "ex : importSNPs_casava('/u/daoudat/py/pyGeno/importationPackages/genomes/ARN_R/ARN_R.tar.gz')"

