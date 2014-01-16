import os, glob, gzip, tarfile, shutil, time
from ConfigParser import SafeConfigParser

import configuration as conf

import rabaDB.setup
#rabaDB.setup._DEBUG_MODE = True
from rabaDB.Raba import *
import rabaDB.fields as rf
from rabaDB.filters import RabaQuery

from Genome import Genome
from Chromosome import Chromosome
from Gene import Gene
from Transcript import Transcript
from Exon import Exon
from Protein import Protein
from SNP import *

from pyGeno.tools.GTFTools import GTFFile

def backUpDB() :
	"backup the current database version. automatically called by importGenome(). Returns the filename of the backup"
	st = time.ctime().replace(' ', '_')
	fn = conf.pyGeno_RABA_DBFILE.replace('.db', '_%s_auto-bck.db' % st)
	shutil.copy2(conf.pyGeno_RABA_DBFILE, fn)

	return fn

def deleteBackUps(forceDelete = False) :
	"if forceDelete = False a confirmation is asked for each file"
	if forceDelete :
		inp = raw_input("WARNING! All backups will be permanently lost, procced? (Y)\n")
		if inp.upper() != 'Y' :
			return

	path = os.path.normpath(conf.pyGeno_SETTINGS_PATH + '/_*_auto-bck.db')
	for f in glob.glob(path) :
		if forceDelete :
			print 'deleting %s...' % f
			os.remove(f)
		else :
			inp = raw_input("delete file %s ? (Y)\n" % f)
			if inp.upper() == 'Y' :
				os.remove(f)
				print '\tdeleted.'

def importGenome(packageFile, verbose = False) :
	r"""Import a pyGeno genome package. A genome packages is a tar.gz ball that contains at it's root:
	-gziped fasta files for all chromosomes
	-gziped GTF gene_set file from ensembl
	-a manifest.ini file such as:
		[package_infos]
		description = Test package. This package installs only chromosome Y of mus musculus
		maintainer = Tariq Daouda
		maintainer_contact = tariq.daouda@umontreal.ca
		version = GRCm38.73

		[genome]
		specie = Mus_musculus
		name = GRCm38_test
		source = http://useast.ensembl.org/info/data/ftp/index.html

		[chromosome_files]
		Y = Mus_musculus.GRCm38.73.dna.chromosome.Y.fa.gz

		[gene_set]
		gtf = Mus_musculus.GRCm38.73_Y-only.gtf.gz

		All files except the manifest can be downloaded at: http://useast.ensembl.org/info/data/ftp/index.html
		A rollback is performed, if an exception is caught during importation
	"""

	def reformatItems(items) :
		s = str(items)
		s = s.replace('[', '').replace(']', '').replace("',", ': ').replace('), ', '\n').replace("'", '').replace('(', '').replace(')', '')
		return s

	print 'importing genome package %s...' % packageFile
	pFile = tarfile.open(packageFile)
	packageDir = os.path.normpath('./.tmp_genome_import')

 	if os.path.isdir(packageDir) :
		shutil.rmtree(packageDir)
	os.makedirs(packageDir)

	for mem in pFile :
		pFile.extract(mem, packageDir)

	parser = SafeConfigParser()
	parser.read(os.path.normpath(packageDir+'/manifest.ini'))
	packageInfos = parser.items('package_infos')

	genomeName = parser.get('genome', 'name')
	specie = parser.get('genome', 'specie')
	genomeSource = parser.get('genome', 'source')
	gtfFile = parser.get('gene_set', 'gtf')
	chromosomesFiles = dict(parser.items('chromosome_files'))
	chromosomeSet = set(chromosomesFiles.keys())

	print "Importing:\n\t%s\nGenome:\n\t%s\n..."  % (reformatItems(packageInfos).replace('\n', '\n\t'), reformatItems(parser.items('genome')).replace('\n', '\n\t'))
	bckFn = backUpDB()
	print "=====\nIf anything goes wrong, the db has been backuped here: %s\nSimply rename it to: %s\n=====" %(bckFn, conf.pyGeno_RABA_DBFILE)

	try :
		genome = Genome(name = genomeName, specie = specie)
		raise ValueError("There seems to be already a genome (%s, %s), please call deleteGenome() first if you want to reinstall it" % (genomeName, specie))
	except KeyError:
		pass

	genome = Genome()
	genome.set(name = genomeName, specie = specie, source = genomeSource, packageInfos = packageInfos)
	seqTargetDir = os.path.normpath(genome.getSequencePath())

	if os.path.isdir(seqTargetDir) :
		raise ValueError("The directory %s already exists, Please call deleteGenome() first if you want to reinstall" % seqTargetDir)
	os.makedirs(seqTargetDir)

	try :
		_importGenomeObjects(os.path.normpath(packageDir+'/'+gtfFile), chromosomeSet, genome, verbose)
		x1Chro = 0
		for chro in genome.chromosomes :
			print "Importing DNA sequence of chromosome %s..." % chro
			length = __importSequence(chro, os.path.normpath(packageDir+'/'+chromosomesFiles[chro.number.lower()]), seqTargetDir)
			chro.x1 = x1Chro
			chro.x2 = x1Chro+length
			x1Chro = chro.x2

		genome.save()
		shutil.rmtree(packageDir)
	except (KeyboardInterrupt, SystemExit, Exception) :
		print "===>Exception caught! Rollback!<==="
		os.remove(conf.pyGeno_RABA_DBFILE)
		os.rename(bckFn, conf.pyGeno_RABA_DBFILE)
		shutil.rmtree(seqTargetDir)
		shutil.rmtree(packageDir)
		raise

def deleteGenome(name, specie) :
	"removes all infos about a genome"
	print 'deleting genome (%s, %s)...' % (specie, name)

	warnings = False
	try :
		shutil.rmtree(conf.getGenomeSequencePath(specie, name))
	except OSError as e:
		print 'WARNING, Unable to delete folder =>', e
		warnings = True

	try :
		print '\tdeleting genome information (%s, %s)...' % (name, specie)
		genome = Genome(name = name, specie = specie)
		for typ in (Chromosome, Gene, Transcript, Exon, Protein) :
			print '\tdeleting all %ss ...' % (typ.__name__)
			try :
				rq = RabaQuery(typ)
				rq.addFilter(genome = genome)
				for e in rq.iterRun() :
					e.delete()
			except OSError as e:
				warnings = True
				print 'WARNING, Unable to delete %s' % typ.__name__

		genome.delete()
	except Exception as e:
		print 'WARNING, Unable to remove genome from db =>', e
		warnings = True

	if warnings :
		print """The deletion could not be fully performed this could be due to various reasons.
		If the problem persist you can try to revert to a previous database version by renaming a pyGenoRaba_[date]_auto-bck.db
		to pyGenoRaba.db in %s and manually erasing folder %s and reinstalling the package.
		""" % (conf.pyGeno_SETTINGS['DATA_PATH'], conf.getGenomeSequencePath(specie, name))

def _importGenomeObjects(gtfFilePath, chroSet, genome, verbose = False) :
	print 'Importing gene set infos from %s...' % gtfFilePath

	gtf = GTFFile()
	gtf.parseFile(gtfFilePath, gziped = True)

	chromosomes = {}
	genes = {}
	transcripts = {}
	proteins = {}
	exons = {}
	for i in range(len(gtf)) :
		chroNumber = str(gtf.get(i, 'seqname'))
		if chroNumber.upper() in chroSet or chroNumber.lower() in chroSet:

			chroNumber = chroNumber.upper()
			if chroNumber not in chromosomes :
				print 'Chromosome %s...' % chroNumber
				chromosomes[chroNumber] = Chromosome()
				chromosomes[chroNumber].set(genome = genome, number = chroNumber)
				chromosomes[chroNumber].dataType = 'heavy'
			try :
				geneId = gtf.get(i, 'gene_id')
				geneName = gtf.get(i, 'gene_name')
			except KeyError :
				if verbose :
					print 'Warning: no gene_id/name found in line %d' % i

			strand = gtf.get(i, 'strand')
			gene_biotype = gtf.get(i, 'gene_biotype')
			regionType = gtf.get(i, 'feature')

			x1 = int(gtf.get(i, 'start')) - 1
			x2 = int(gtf.get(i, 'end'))
			if x1 > x2 :
				x1, x2 = x2, x1

			if geneId not in genes :
				#if verbose :
				print '\tGene %s, %s...' % (geneId, geneName)
				genes[geneId] = Gene()
				genes[geneId].set(genome = genome, id = geneId, chromosome = chromosomes[chroNumber], name = geneName, strand = strand, biotype = gene_biotype)

			try :
				transId = gtf.get(i, 'transcript_id')
				transName = gtf.get(i, 'transcript_name')
			except KeyError :
				if verbose :
					print '\tWarning: no transcript_id, name found in line %d' % i

			if transId not in transcripts :
				if verbose :
					print '\tTranscript %s, %s...' % (transId, transName)
				transcripts[transId] = Transcript(importing = True)
				transcripts[transId].set(genome = genome, id = transId, chromosome = chromosomes[chroNumber], gene = genes[geneId], name = transName)
			try :
				protId = gtf.get(i, 'protein_id')
				if protId not in proteins :
					if verbose :
						print '\tProtein %s...' % (protId)
					proteins[protId] = Protein(importing = True)
					proteins[protId].set(genome = genome, id = protId, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], name = transName)
					transcripts[transId].protein = proteins[protId]
					#proteins[protId].save()
			except KeyError :
				if verbose :
					print 'Warning: no protein_id found in line %d' % i

			try :
				exonNumber = gtf.get(i, 'exon_number')
			except KeyError :
				print 'Warning: no number found in line %d' % i

			exonKey = (transId, exonNumber)

			if regionType == 'exon' :
				try :
					exonId = gtf.get(i, 'exon_id')
					if exonId not in exons :
						exons[exonKey] = Exon(importing = True)
						exons[exonKey].set(genome = genome, id = exonId, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], strand = strand, number = exonNumber, x1 = x1, x2 = x2)
						#exons[exonKey].save()
				except KeyError :
					print 'Warning: no exon_id found in line %d' % i

			elif regionType == 'CDS' :
				exons[exonKey].CDS_x1 = x1
				exons[exonKey].CDS_x2 = x2
			#elif regionType == 'start_codon' :
			#	exons[exonKey].startCodon = x1
			elif regionType == 'stop_codon' :
				if strand == '+' :
					exons[exonKey].x2 += 3
					if exons[exonKey].CDS_x2 != None :
						exons[exonKey].CDS_x2 += 3
				if strand == '-' :
					exons[exonKey].x1 -= 3
					if exons[exonKey].CDS_x1 != None :
						exons[exonKey].CDS_x1 -= 3

	print 'creating relations...'
	print '\ttranscript.exons...'
	for transcript in transcripts.itervalues() :
		f = RabaQuery(Exon)
		f.addFilter(**{'transcript' : transcript})
		p = f.run()
		transcript.exons = f.run()

	print '\tgene.transcripts/.exons...'
	for gene in genes.itervalues() :
		f = RabaQuery(Transcript)
		f.addFilter(**{'gene' : gene})
		gene.transcripts = f.run()

		f = RabaQuery(Exon)
		f.addFilter(**{'gene' : gene})
		gene.exons = f.run()

	print '\tchromosome.genes...'
	for chro in chromosomes.itervalues() :
		f = RabaQuery(Gene)
		f.addFilter(**{'chromosome' : chro})
		chro.genes = f.run()

	print '\tgenome.chromosomes...'
	genome.chromosomes = RabaList(chromosomes.values())
	print 'saving...'
	print 'Done.'
	genome.save()
	#print "=========", transcript.exons
	try :
		refGenomeName = conf.getReferenceGenome(genome.specie)
	except KeyError:
		refGenomeName = genome.name
		conf.setReferenceGenome(genome.specie, genome.name)

	print 'Current reference genome for specie %s is %s' %(genome.specie, genome.name)

def __importSequence(chromosome, fastaFile, targetDir) :
	print 'making data for chromsome %s, source file: %s...' %(chromosome.number, fastaFile)

	f = gzip.open(fastaFile)
	header = f.readline()
	strRes = f.read().upper().replace('\n', '').replace('\r', '')
	f.close()

	fn = '%s/chromosome%s.dat' % (targetDir, chromosome.number)
	f = open(fn, 'w')
	f.write(strRes)
	f.close()
	chromosome.dataFile = fn
	chromosome.header = header
	return len(strRes)

def importGenome_casava(specie, genomeName, snpsTxtFile) :
	"""TODO: TEST"""

	print 'importing %s genome %s...' % (specie, genomeName)

	genome = Genome(name = genomeName, specie = specie)
	genome.name = genomeName
	genome.specie = specie
	genome.source = snpsTxtFile

	f = open(snpsTxtFile)
	lines = f.readlines()
	f.close()

	currChrNumber = '-1'
	for l in lines :
		if l[0] != '#' : #ignore comments
			sl = l.replace('\t\t', '\t').split('\t')
			if sl[0] != currChrNumber :
				print 'importing snp data for chromosome %s...' % sl[0]
				currChrNumber = sl[0]
				chromosome = Chromosome()
				chromosome.number = currChrNumber
				chromosome.genome = genome
				chromosome.dataType = 'CasavaSNP'

			snp = Casava_SNP()
			snp.chromosome = chromosome
			snp.genome = genome
			#first column: chro, second first of range (identical to second column)
			snp.values['pos'] = int(sl[2])
			snp.values['bcalls_used'] = sl[3]
			snp.values['bcalls_filt'] = sl[4]
			snp.values['ref'] = sl[5]
			snp.values['QSNP'] = int(sl[6])
			snp.values['max_gt'] = uf.getPolymorphicNucleotide(sl[7])
			snp.values['Qmax_gt'] = int(sl[8])
			snp.values['max_gt_poly_site'] = sl[9]
			snp.values['Qmax_gt_poly_site'] = int(sl[10])
			snp.values['A_used'] = int(sl[11])
			snp.values['C_used'] = int(sl[12])
			snp.values['G_used'] = int(sl[13])
			snp.values['T_used'] = int(sl[14])

	print 'saving...'
	genome.save()
	print 'importation %s of genome %s done.' %(specie, genomeName)


def import_dbSNP(packageFolder, specie, versionName) :
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

	RabaConnection(conf.pyGeno_RABA_NAMESPACE).autoOnSaveCommit(False)

	def parseSNP(snpLines, specie, chroNumber, version) :
		lines = snpLines.split('\n')

		snp = dbSNP_SNP()
		snp.version = version
		snp.specie = specie
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

		#snp.save()

	files = glob.glob(os.path.normpath(packageFolder+'/*.flat.gz'))

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
			parseSNP(snp, specie, chroNumber, versionName)

	print 'saving...'
	RabaConnection(conf.pyGeno_RABA_NAMESPACE).commit()
	RabaConnection(conf.pyGeno_RABA_NAMESPACE).autoOnSaveCommit(True)
	print 'done.'

if __name__ == "__main__" :
	#deleteBackUps(forceDelete = True)
	deleteGenome(specie = 'human', name = 'GRCh37.74')
	#importGenome('/u/daoudat/py/pyGeno/importationPackages/genomes/mouse/mus_musculus_Y-only.tar.gz', verbose = False)
	importGenome('/u/daoudat/py/pyGeno/importationPackages/GRCh37.74/GRCh37.74.tar.gz', verbose = False)

	#g = Genome(specie = 'Mus_musculus', name = 'GRCm38_test')
	#a = g.get(Gene)
	#print a[0].transcripts#, a[0].help()
	#a[0].genome.chromosomes = [a[0]]
	#g.save()

	#importGenome_casava('Mus_musculus', 'musCasava', 'http://www.bioinfo.iric.ca/seq/Project_DSP008a/Build_Diana_ARN_R/snps.txt')

	#import_dbSNP('/u/daoudat/py/pyGeno/importationPackages/dbSNP/human/test', 'Mus_musculus', 'test')
	#print g.chromosomes
	#f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, dbSNP_SNP)
	#f.addFilter(**{'version' : 'test'})
	#for e in f.run() :
	#	print e.rsId
