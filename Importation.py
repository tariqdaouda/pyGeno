import os, glob, gzip, tarfile, shutil, time, sys, gc, cPickle, tempfile
from ConfigParser import SafeConfigParser

import configuration as conf

import rabaDB.setup
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

import cProfile

def printf(*s) :
	for e in s[:-1] :
		print e,
	print s[-1]

	sys.stdout.flush()

def backUpDB() :
	"backup the current database version. automatically called by importGenome(). Returns the filename of the backup"
	st = time.ctime().replace(' ', '_')
	fn = conf.pyGeno_RABA_DBFILE.replace('.db', '_%s_auto-bck.db' % st)
	shutil.copy2(conf.pyGeno_RABA_DBFILE, fn)

	return fn

"""
def deleteBackUps(forceDelete = False) :
	"if forceDelete = False a confirmation is asked for each file"
	if forceDelete :
		inp = raw_input("WARNING! All backups will be permanently lost, procced? (Y)\n")
		if inp.upper() != 'Y' :
			return

	path = os.path.normpath(conf.pyGeno_SETTINGS_PATH + '/_*_auto-bck.db')
	for f in glob.glob(path) :
		if forceDelete :
			printf('deleting %s...' % f)
			os.remove(f)
		else :
			inp = raw_input("delete file %s ? (Y)\n" % f)
			if inp.upper() == 'Y' :
				os.remove(f)
				printf('\tdeleted.')
"""
def _decompressPackage(packageFile) :
	pFile = tarfile.open(packageFile)
	packageDir = os.path.normpath('%s/tmp_pygeno_import' % tempfile.gettempdir())

 	if os.path.isdir(packageDir) :
		shutil.rmtree(packageDir)
	os.makedirs(packageDir)

	for mem in pFile :
		pFile.extract(mem, packageDir)

	return packageDir

#def _checkGenomeNotInDb(genomeName, specie) :
#	"raise an excpetion is the genome is already present in the db"
#	try :
#		genome = Genome(name = genomeName, specie = specie)
#		raise ValueError("There seems to be already a genome (%s, %s), please call deleteGenome() first if you want to reinstall it" % (genomeName, specie))
#	except KeyError:
#		pass

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

	printf('importing genome package %s...' % packageFile)

	packageDir = _decompressPackage(packageFile)

	parser = SafeConfigParser()
	parser.read(os.path.normpath(packageDir+'/manifest.ini'))
	packageInfos = parser.items('package_infos')

	genomeName = parser.get('genome', 'name')
	specie = parser.get('genome', 'specie')
	genomeSource = parser.get('genome', 'source')
	gtfFile = parser.get('gene_set', 'gtf')
	chromosomesFiles = dict(parser.items('chromosome_files'))
	chromosomeSet = set(chromosomesFiles.keys())

	try :
		genome = Genome(name = genomeName, specie = specie)
		#deleteGenome(specie, genomeName)
		raise ValueError("There seems to be already a genome (%s, %s), please call deleteGenome() first if you want to reinstall it" % (genomeName, specie))
	except KeyError:
		pass

	seqTargetDir = conf.getGenomeSequencePath(specie, genomeName)
	if os.path.isdir(seqTargetDir) :
		raise ValueError("The directory %s already exists, Please call deleteGenome() first if you want to reinstall" % seqTargetDir)

	#try :
	os.makedirs(seqTargetDir)
	#except OSError:
	#	printf("Warning the directory %s already exists" % seqTargetDir)

	genome = Genome()
	genome.set(name = genomeName, specie = specie, source = genomeSource, packageInfos = packageInfos)

	#try :
	printf("Importing:\n\t%s\nGenome:\n\t%s\n..."  % (reformatItems(packageInfos).replace('\n', '\n\t'), reformatItems(parser.items('genome')).replace('\n', '\n\t')))
	#bckFn = backUpDB()
	#printf("=====\nIf anything goes wrong, the db has been backuped here: %s\nSimply rename it to: %s\n=====" %(bckFn, conf.pyGeno_RABA_DBFILE))

	chros = _importGenomeObjects(os.path.normpath(packageDir+'/'+gtfFile), chromosomeSet, genome, verbose)
	x1Chro = 0
	for chro in chros :
		printf("Importing DNA sequence of chromosome %s..." % chro)
		length = _importSequence(chro, os.path.normpath(packageDir+'/'+chromosomesFiles[chro.number.lower()]), seqTargetDir)
		chro.x1 = x1Chro
		chro.x2 = x1Chro+length
		x1Chro = chro.x2

	genome.save()
	shutil.rmtree(packageDir)
	#except (KeyboardInterrupt, SystemExit, Exception) :
	#	raise
	#	printf("===>Exception caught! Rollback!<===")
	#	os.remove(conf.pyGeno_RABA_DBFILE)
	#	os.rename(bckFn, conf.pyGeno_RABA_DBFILE)
	#	shutil.rmtree(seqTargetDir)
	#	shutil.rmtree(packageDir)
	#	raise

def deleteGenome(specie, name) :
	"removes all infos about a genome"

	printf('deleting genome (%s, %s)...' % (specie, name))

	rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).beginTransaction()
	printf('\tdeleting genome information (%s, %s)...' % (specie, name))
	objs = []
	genome = Genome(name = name, specie = specie)
	objs.append(genome)
	for typ in (Chromosome, Gene, Transcript, Exon, Protein) :
		printf('\tgathering %ss for deletion...' % (typ.__name__))
		for e in genome.get(typ) :
			objs.append(e)

	printf('\tdeleting...')
	for e in objs :
		e.delete()

	try :
		shutil.rmtree(conf.getGenomeSequencePath(specie, name))
	except OSError as e:
		printf('WARNING, Unable to delete folder: ', e)

	rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).endTransaction()

#@profile
def _importGenomeObjects(gtfFilePath, chroSet, genome, verbose = 0) :
	"verbose is int [0, 4] for various levels of verbosity"

	printf('Importing gene set infos from %s...' % gtfFilePath)
	startTime = time.time()

	gtf = GTFFile()
	gtf.parseFile(gtfFilePath, gziped = True)

	chromosomes = {}
	genes = {}
	transcripts = {}
	proteins = {}
	exons = {}
	chroNumber = None
	chroStartTime = time.time()
	for i in range(len(gtf)) :
		chroN = str(gtf.get(i, 'seqname'))

		if (chroN.upper() in chroSet or chroN.lower() in chroSet):
			if chroN != chroNumber and chroNumber != None :
				printf('\tdone (%fmin), total time (%f).' %((time.time()-chroStartTime)/60, (time.time()-startTime)/60))
				chroStartTime = time.time()

			chroNumber = chroN.upper()
			if chroNumber not in chromosomes :
				printf('Importing objects for chromosome %s...' % chroNumber)
				chromosomes[chroNumber] = Chromosome(importing = True)
				chromosomes[chroNumber].set(genome = genome, number = chroNumber)
				chromosomes[chroNumber].dataType = 'heavy'
			try :
				geneId = gtf.get(i, 'gene_id')
				geneName = gtf.get(i, 'gene_name')
			except KeyError :
				if verbose :
					printf('Warning: no gene_id/name found in line %d' % i)

			strand = gtf.get(i, 'strand')
			gene_biotype = gtf.get(i, 'gene_biotype')
			regionType = gtf.get(i, 'feature')

			x1 = int(gtf.get(i, 'start')) - 1
			x2 = int(gtf.get(i, 'end'))
			if x1 > x2 :
				x1, x2 = x2, x1

			if geneId not in genes :
				if verbose > 0 :
					printf('\tGene %s, %s...' % (geneId, geneName))
				genes[geneId] = Gene()
				genes[geneId].set(genome = genome, id = geneId, chromosome = chromosomes[chroNumber], name = geneName, strand = strand, biotype = gene_biotype)

			try :
				transId = gtf.get(i, 'transcript_id')
				transName = gtf.get(i, 'transcript_name')
			except KeyError :
				if verbose > 1 :
					printf('\t\tWarning: no transcript_id, name found in line %d' % i)

			if transId not in transcripts :
				if verbose > 1 :
					printf('\t\tTranscript %s, %s...' % (transId, transName))
				transcripts[transId] = Transcript(importing = True)
				transcripts[transId].set(genome = genome, id = transId, chromosome = chromosomes[chroNumber], gene = genes[geneId], name = transName)
			try :
				protId = gtf.get(i, 'protein_id')
				if protId not in proteins :
					if verbose > 1 :
						printf('\t\tProtein %s...' % (protId))
					proteins[protId] = Protein(importing = True)
					proteins[protId].set(genome = genome, id = protId, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], name = transName)
					transcripts[transId].protein = proteins[protId]

			except KeyError :
				if verbose > 2 :
					printf('Warning: no protein_id found in line %d' % i)

			try :
				exonNumber = gtf.get(i, 'exon_number')
			except KeyError :
				printf('Warning: no number found in line %d' % i)

			exonKey = (transId, exonNumber)

			if regionType == 'exon' :
				try :
					exonId = gtf.get(i, 'exon_id')
					if verbose > 3 :
						printf('\t\t\texon %s...' % (exonId))
					if exonId not in exons :
						exons[exonKey] = Exon(importing = True)
						exons[exonKey].set(genome = genome, id = exonId, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], strand = strand, number = exonNumber, x1 = x1, x2 = x2)
						transcripts[transId].exons.append(exons[exonKey])

				except KeyError :
					printf('Warning: no exon_id found in line %d' % i)

			elif regionType == 'CDS' :
				exons[exonKey].CDS_x1 = x1
				exons[exonKey].CDS_x2 = x2
			elif regionType == 'stop_codon' :
				if strand == '+' :
					exons[exonKey].x2 += 3
					if exons[exonKey].CDS_x2 != None :
						exons[exonKey].CDS_x2 += 3
				if strand == '-' :
					exons[exonKey].x1 -= 3
					if exons[exonKey].CDS_x1 != None :
						exons[exonKey].CDS_x1 -= 3

	printf('\tdone (%fmin), total time (%f).' %((time.time()-chroStartTime)/60, (time.time()-startTime)/60))

	rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).beginTransaction()
	if verbose > 0 :
		printf('saving %d chromsomes...' % len(chromosomes))
	for c in chromosomes.itervalues() :
		c.save()

	if verbose > 0 :
		printf('saving %d genes...' % len(genes))
	for c in genes.itervalues() :
		c.save()

	if verbose > 0 :
		printf('saving %d transcripts and exons %d...' % (len(transcripts), len(exons)))
	for c in transcripts.itervalues() :
		c.save()

	if verbose > 0 :
		printf('saving %d proteins...' % len(proteins))
	for c in proteins.itervalues() :
		c.save()
	genome.save()
	printf('\tcommiting changes...')
	rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).endTransaction()
	printf('\tdone total time (%f).' %((time.time()-startTime)/60))

	try :
		refGenomeName = conf.getReferenceGenome(genome.specie)
	except KeyError:
		refGenomeName = genome.name
		conf.setReferenceGenome(genome.specie, genome.name)
		printf('Auto setting: Current reference genome for specie %s is now %s' %(genome.specie, genome.name))

	return chromosomes.values()

def _importSequence(chromosome, fastaFile, targetDir) :
	printf('making data for chromsome %s, source file: %s...' %(chromosome.number, fastaFile))

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

def importSNPs_casava(specie, setName, packageFile) :
	"""TODO: TEST"""

	printf('importing %s genome %s...' % (specie, genomeName))
	#bckFn = backUpDB()
	printf("=====\nIf anything goes wrong, the db has been backuped here: %s\nSimply rename it to: %s\n=====" %(bckFn, conf.pyGeno_RABA_DBFILE))

	packageDir = _decompressPackage(packageFile)

	parser = SafeConfigParser()
	parser.read(os.path.normpath(packageDir+'/manifest.ini'))
	packageInfos = parser.items('package_infos')

	genomeName = parser.get('genome', 'name')
	specie = parser.get('genome', 'specie')
	genomeSource = parser.get('genome', 'source')
	snpsTxtFile = os.path.normpath('%s/%s' %(packageDir, parser.get('snps', 'casava_snps_txt')))

	#try :
	rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).beginTransaction()

	f = open(snpsTxtFile)
	lines = f.readlines()
	f.close()

	currChrNumber = None
	for l in lines :
		if l[0] != '#' : #ignore comments
			sl = l.replace('\t\t', '\t').split('\t')
			if sl[0] != currChrNumber :
				printf('importing snp data for chromosome %s...' % sl[0])
				currChrNumber = sl[0].replace('CHR', '').upper()

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

	snpMaster = SNPMaster()
	snpMaster.set(setName = setName, SNPType = 'CasavaSNP', specie = specie)
	snpMaster.save()

	printf('saving...')
	rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).endTransaction()
	printf('creating indexes...')
	CasavaSNP.addIndex('setName')
	CasavaSNP.addIndex('start')
	printf('importation %s of genome %s done.' %(specie, genomeName))

	#except (KeyboardInterrupt, SystemExit, Exception) :
	#	printf("===>Exception caught! Rollback!<===")
	#	os.remove(conf.pyGeno_RABA_DBFILE)
	#	os.rename(bckFn, conf.pyGeno_RABA_DBFILE)
	#	raise

if __name__ == "__main__" :
	#deleteBackUps(forceDelete = True)(, )
	#importGenome('/u/daoudat/py/pyGeno/importationPackages/genomes/mouse/mus_musculus_Y-only.tar.gz', verbose = 0)
	#deleteGenome(specie = 'Mus_musculus', name = 'GRCm38_test')
	#deleteGenome(specie = 'human', name = 'GRCh37.74')
	importGenome('/u/daoudat/py/pyGeno/importationPackages/GRCh37.74/GRCh37.74.tar.gz', verbose = 0)
	#g = Genome(specie = 'Mus_musculus', name = 'GRCm38_test')
	#printf(g.get(Transcript)[0].exons)
	#deleteGenome(specie = 'human', name = 'ARN_R')
	#importGenome_casava('human', 'ARN_R', '/u/daoudat/py/pyGeno/importationPackages/genomes/ARN_R/ARN_R.tar.gz')
	#importGenome_casava('human', 'ARN_R', '/u/corona/Project_DSP008a/Build_Diana_ARN_R/snps.txt')
	#g = Genome(specie = 'human', name = 'ARN_R')
	#printf(g.snps)
	#printf(Genome(specie = 'human', name = 'GRCh37.74'))

