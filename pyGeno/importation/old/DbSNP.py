#dbSNP ASN.1 parser, may be broken

import configuration as conf

from rabaDB.setup import *
RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
from rabaDB.Raba import *
import rabaDB.fields as rf

from tools import UsefulFunctions as uf
from exceptions import *

def import_dbSNP(packageFile) :
	"""To import dbSNP informations, download ASN1_flat files from the
	dbSNP ftp : ftp://ftp.ncbi.nih.gov/snp/organisms/ and place them all in one single folder. This folder
	will be considered as a package.
	Launch this function and go make yourself a cup of coffee, this function has absolutly not been written to be fast

	versionName is name with wich you want to call this specific version of dbSNP

	pyGeno snp format is somewhat different from dbSNP's :
		- all snps have an orientation of +. If a snp has an orientation of -, it's alleles are replaced by their complements
		- positions are 0 based
		- the only value extracted from the files are : 'posistion', 'rs', 'type', 'assembly', 'chromosome', 'validated', 'alleles', 'original_orientation', 'maf_allele', 'maf_count', 'maf', 'het', 'se(het)
		-loc is a dictionary allele wise that simplifies the line loc'
	"""
	#TODO: make it a tar ball package with a manifest.ini

	def parseSNP(snpLines, chroNumber) :
		lines = snpLines.split('\n')

		snp = dbSNP_SNP()
		snp.chromosomeNumber = chroNumber

		#numericFields: maf_count', 'maf', 'het', 'se(het)'
		for l in lines :
			sl = l.split('|')
			if sl[0][:2] == 'rs' :
				snp.rsId = sl[0][2:].strip()
				snp.type = sl[3].strip()

			elif sl[0][:3] == 'SNP' :
				snp.alleles = sl[1].strip().replace('alleles=', '').replace("'", "")
				het = sl[2].strip().replace('het=', '')
				try :
					snp.het = float(het)
				except :
					pass

				se_het = sl[3].strip().replace('se(het)=', '')
				try :
					snp.se_het = float(se_het)
				except :
					pass

			elif sl[0][:3] == 'VAL' :
				snp.validated = sl[1].strip().replace("validated=", '')

			elif sl[0][:3] == 'CTG' and sl[1].find('GRCh') > -1 :
				snp.original_orientation = sl[-1].replace('orient=', '').strip()
				snp.assembly = sl[1].replace('assembly=', '').strip()
				snp.chromosome = sl[2].replace('chr=', '').strip()
				pos = sl[3].replace('chr-pos=', '').strip()

				try:
					snp.pos = int(pos) -1
				except :
					snp.pos = pos

			elif sl[0][:4] == 'GMAF' :
				snp.maf_allele = sl[1].strip().replace('allele=', '')
				maf_count = sl[2].strip().replace('count=', '')
				try :
					snp.maf_count = float(maf_count)
				except :
					snp.maf_count = maf_count

				maf = sl[3].strip().replace('MAF=', '')
				try :
					snp.maf = float(maf)
				except :
					snp.maf = maf

			elif sl[0][:3] == 'LOC' :
				loc = dbSNP_SNPLOC()
				loc.allele = sl[4].strip().replace('allele=', '')
				loc.gene = sl[1].strip()
				loc.fxn_class = sl[3].strip().replace('fxn-class=', '')
				try :
					loc.residue = sl[6].strip().replace('residue=', '')
				except IndexError :
					pass

		if snp.original_orientation == '-' :
			snp.alleles = uf.complement(snp.alleles)

		return snp

	try :
		dbSNP = DbSNP(specie = specie, version = version)
		raise ValueError("There seems to be already a genome (%s, %s), please call deleteDbSNP() first if you want to reinstall it" % (genomeName, specie))
	except KeyError:
		pass

	dbSNP = DbSNP()
	dbSNP.set(specie = specie, version = version)

	print 'importing dbSNP package %s...' % packageFile
	pFile = tarfile.open(packageFile)
	packageDir = os.path.normpath('./.tmp_dbsnp_import')

	if os.path.isdir(packageDir) :
		shutil.rmtree(packageDir)
	os.makedirs(packageDir)

	for mem in pFile :
		pFile.extract(mem, packageDir)

	parser = SafeConfigParser()
	parser.read(os.path.normpath(packageDir+'/manifest.ini'))
	packageInfos = parser.items('package_infos')

	genomeName = parser.get('dbSNP', 'specie')
	specie = parser.get('dbSNP', 'version')
	genomeSource = parser.get('dbSNP', 'source')
	gtfFile = parser.get('gene_set', 'gtf')
	chromosomesFiles = dict(parser.items('chromosome_files'))
	chromosomeSet = set(chromosomesFiles.keys())

	print "Importing:\n\t%s\nGenome:\n\t%s\n..."  % (reformatItems(packageInfos).replace('\n', '\n\t'), reformatItems(parser.items('genome')).replace('\n', '\n\t'))
	bckFn = backUpDB()
	print "=====\nIf anything goes wrong, the db has been backuped here: %s\nSimply rename it to: %s\n=====" %(bckFn, conf.pyGeno_RABA_DBFILE)

	files = glob.glob(os.path.normpath(packageDir+'/*.flat.gz'))

	for fil in files :
		chrStrStartPos = fil.find('ch')
		chrStrStopPos = fil.find('.flat')

		chroNumber = fil[chrStrStartPos+2: chrStrStopPos]

		print "extracting file :", fil, "..."
		f = gzip.open(fil)
		snps = f.read().split('\n\n')
		f.close()

		print "\timporting snps..."
		for snp in snps[1:] :
			dbSNP.snps.append(parseSNP(snp, specie, chroNumber, versionName))

	print 'saving...'
	dbSNP.save()
	print 'done.'

class DbSNP(Raba) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	version = rf.Primitive()
	snps = rf.Relation('dbSNP_SNP')

	def __init__(self) :
		pass

class dbSNP_SNP(pyGenoObject) :

	rsId = rf.Primitive()
	pos = rf.Primitive()
	type = rf.Primitive()
	alleles = rf.Primitive()
	validated = rf.Primitive()

	assembly = rf.Primitive()
	original_orientation = rf.Primitive()
	maf_allele = rf.Primitive()
	maf_count = rf.Primitive()
	maf = rf.Primitive()
	het = rf.Primitive()
	se_het = rf.Primitive()

	loc = rf.RList()

	specie = rf.Primitive()
	chromosomeNumber = rf.Primitive()

	def __init__(self, *args, **fieldsSet) :
		pass

	def _curate(self) :
		pass

class dbSNP_SNPLOC(pyGenoObject) :

	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	allele = rf.Primitive()
	fxn_class = rf.Primitive()
	gene = rf.Primitive()
	residue = rf.Primitive()

	snp = rf.RabaObject('dbSNP_SNP')

	def __init__(self) :
		pass

	def _curate(self) :
		pass
