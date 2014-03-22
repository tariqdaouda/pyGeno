import os, glob, gzip, tarfile, shutil, time, sys, gc, cPickle, tempfile
from ConfigParser import SafeConfigParser

import configuration as conf

from SNP import *

from pyGeno.tools.ProgressBar import ProgressBar
from pyGeno.tools.io import printf
from Genomes import _decompressPackage


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

	f = open(snpsTxtFile)
	lines = f.readlines()
	f.close()

	CasavaSNP.dropIndex('setName')
	CasavaSNP.dropIndex('start')

	pBar = ProgressBar(len(lines))
	pLabel = ''
	currChrNumber = None
	for l in lines :
		if l[0] != '#' : #ignore comments
			sl = l.replace('\t\t', '\t').split('\t')
			tmpChr = sl[0].upper().replace('CHR', '')
			if tmpChr != currChrNumber :
				currChrNumber = tmpChr
				pLabel = 'Chr %s...' % currChrNumber

			snp = CasavaSNP()
			snp.chromosomeNumber = currChrNumber
			snp.specie = specie
			snp.setName = setName
			#first column: chro, second first of range (identical to second column)
			snp.start = int(sl[2])
			snp.end = int(sl[2])+1
			snp.bcalls_used = sl[3]
			snp.bcalls_filt = sl[4]
			snp.ref = sl[5]
			snp.QSNP = int(sl[6])
			snp.alleles = uf.getPolymorphicNucleotide(sl[7]) #max_gt
			snp.Qmax_gt = int(sl[8])
			snp.max_gt_poly_site = sl[9]
			snp.Qmax_gt_poly_site = int(sl[10])
			snp.A_used = int(sl[11])
			snp.C_used = int(sl[12])
			snp.G_used = int(sl[13])
			snp.T_used = int(sl[14])
			snp.save()
			pBar.update(label = pLabel)

	snpMaster = SNPMaster()
	snpMaster.set(setName = setName, SNPType = 'CasavaSNP', specie = specie)
	snpMaster.save()

	printf('saving...')
	conf.db.endTransaction()
	printf('creating indexes...')
	CasavaSNP.requireIndex('setName')
	CasavaSNP.requireIndex('start')
	printf('importation of SNP set %s for specie %s done.' %(setName, specie))

if __name__ == "__main__" :
	print "ex : importSNPs_casava('/u/daoudat/py/pyGeno/importationPackages/genomes/ARN_R/ARN_R.tar.gz')"

