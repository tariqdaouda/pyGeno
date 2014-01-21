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

	rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).setAutoCommit(False)

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


	try :
		genome = Genome(name = genomeName, specie = specie)
		raise ValueError("There seems to be already a genome (%s, %s), please call deleteGenome() first if you want to reinstall it" % (genomeName, specie))
	except KeyError:
		pass

	seqTargetDir = conf.getGenomeSequencePath(specie, genomeName)
	if os.path.isdir(seqTargetDir) :
		raise ValueError("The directory %s already exists, Please call deleteGenome() first if you want to reinstall" % seqTargetDir)
	os.makedirs(seqTargetDir)

	genome = Genome()
	genome.set(name = genomeName, specie = specie, source = genomeSource, packageInfos = packageInfos)


	try :
		print "Importing:\n\t%s\nGenome:\n\t%s\n..."  % (reformatItems(packageInfos).replace('\n', '\n\t'), reformatItems(parser.items('genome')).replace('\n', '\n\t'))
		bckFn = backUpDB()
		print "=====\nIf anything goes wrong, the db has been backuped here: %s\nSimply rename it to: %s\n=====" %(bckFn, conf.pyGeno_RABA_DBFILE)

		_importGenomeObjects(os.path.normpath(packageDir+'/'+gtfFile), chromosomeSet, genome, verbose)
		x1Chro = 0
		for chro in genome.get(Chromosome) :
			print "Importing DNA sequence of chromosome %s..." % chro
			length = _importSequence(chro, os.path.normpath(packageDir+'/'+chromosomesFiles[chro.number.lower()]), seqTargetDir)
			chro.x1 = x1Chro
			chro.x2 = x1Chro+length
			x1Chro = chro.x2

		genome.save()
		rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).commit()
		rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).setAutoCommit(True)
		shutil.rmtree(packageDir)
	except (KeyboardInterrupt, SystemExit, Exception) :
		print "===>Exception caught! Rollback!<==="
		os.remove(conf.pyGeno_RABA_DBFILE)
		os.rename(bckFn, conf.pyGeno_RABA_DBFILE)
		shutil.rmtree(seqTargetDir)
		shutil.rmtree(packageDir)
		raise

def deleteGenome(specie, name) :
	"removes all infos about a genome"
	print 'deleting genome (%s, %s)...' % (specie, name)

	warnings = False
	rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).setAutoCommit(False)
	try :
		print '\tdeleting genome information (%s, %s)...' % (specie, name)
		genome = Genome(name = name, specie = specie)
		genome.delete()
		for typ in (Chromosome, Gene, Transcript, Exon, Protein) :
			print '\tdeleting all %ss ...' % (typ.__name__)
			try :
				rq = RabaQuery(typ)
				rq.addFilter(genome = genome)
				for e in rq.iterRun() :
					e.delete()
					rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).commit()
			except OSError as e:
				warnings = True
				print 'WARNING, Unable to delete %s' % typ.__name__

	except Exception as e:
		print 'WARNING, Unable to remove Genome from db =>', e
		warnings = True

	try :
		shutil.rmtree(conf.getGenomeSequencePath(specie, name))
	except OSError as e:
		print 'WARNING, Unable to delete folder =>', e
		warnings = True

	rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).commit()
	rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).setAutoCommit(True)

	if warnings :
		print """If the deletion was not fully performed and the problem persist, you can try to revert to a previous database version by renaming a pyGenoRaba_[date]_auto-bck.db
to pyGenoRaba.db in %s and manually erasing folder %s and reinstalling the package.
		""" % (conf.pyGeno_SETTINGS['DATA_PATH'], conf.getGenomeSequencePath(specie, name))

def _importGenomeObjects(gtfFilePath, chroSet, genome, verbose = 0) :
	"verbose is int [0, 4] for various levels of verbosity"

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
				if verbose > 0 :
					print '\tGene %s, %s...' % (geneId, geneName)
				genes[geneId] = Gene()
				genes[geneId].set(genome = genome, id = geneId, chromosome = chromosomes[chroNumber], name = geneName, strand = strand, biotype = gene_biotype)
				#chromosomes[chroNumber].genes.append(genes[geneId])
				#genes[geneId].save()
			try :
				transId = gtf.get(i, 'transcript_id')
				transName = gtf.get(i, 'transcript_name')
			except KeyError :
				if verbose > 1 :
					print '\t\tWarning: no transcript_id, name found in line %d' % i

			if transId not in transcripts :
				if verbose > 1 :
					print '\t\tTranscript %s, %s...' % (transId, transName)
				transcripts[transId] = Transcript(importing = True)
				transcripts[transId].set(genome = genome, id = transId, chromosome = chromosomes[chroNumber], gene = genes[geneId], name = transName)
				#genes[geneId].transcripts.append(transcripts[transId])
				#transcripts[transId].save()
			try :
				protId = gtf.get(i, 'protein_id')
				if protId not in proteins :
					if verbose > 1 :
						print '\t\tProtein %s...' % (protId)
					proteins[protId] = Protein(importing = True)
					proteins[protId].set(genome = genome, id = protId, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], name = transName)
					transcripts[transId].protein = proteins[protId]
					#proteins[protId].save()
			except KeyError :
				if verbose > 2 :
					print 'Warning: no protein_id found in line %d' % i

			try :
				exonNumber = gtf.get(i, 'exon_number')
			except KeyError :
				print 'Warning: no number found in line %d' % i

			exonKey = (transId, exonNumber)

			if regionType == 'exon' :
				try :
					exonId = gtf.get(i, 'exon_id')
					if verbose > 3 :
						print '\t\t\texon %s...' % (exonId)
					if exonId not in exons :
						exons[exonKey] = Exon(importing = True)
						exons[exonKey].set(genome = genome, id = exonId, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], strand = strand, number = exonNumber, x1 = x1, x2 = x2)
						transcripts[transId].exons.append(exons[exonKey])
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
	"""
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
	"""
	print 'saving...'
	"""print '\tgenome...'
	genome.save()
	print '\tchromosomes...'
	for c in chromosomes.itervalues() :
		c.save()
	print '\tgenes...'
	for c in genes.itervalues() :
		c.save()
	print '\ttranscripts + proteins + exons...'"""
	for c in transcripts.itervalues() :
		c.save()
	#print '\tproteins...'
	#for c in proteins.itervalues :
	#	c.save()
	#print '\texons...'
	#for c in exons.itervalues :
	#	c.save()
	#print 'Done.'

	#print "=========", transcript.exons
	try :
		refGenomeName = conf.getReferenceGenome(genome.specie)
	except KeyError:
		refGenomeName = genome.name
		conf.setReferenceGenome(genome.specie, genome.name)

	print 'Current reference genome for specie %s is %s' %(genome.specie, genome.name)

def _importSequence(chromosome, fastaFile, targetDir) :
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
	bckFn = backUpDB()
	print "=====\nIf anything goes wrong, the db has been backuped here: %s\nSimply rename it to: %s\n=====" %(bckFn, conf.pyGeno_RABA_DBFILE)

	try :
		rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).setAutoCommit(False)

		try :
			genome = Genome(name = genomeName, specie = specie)
			raise ValueError("There seems to be already a genome (%s, %s), please call deleteGenome() first if you want to reinstall it" % (genomeName, specie))
		except KeyError:
			pass

		genome = Genome()
		genome.set(name = genomeName, specie = specie, source = snpsTxtFile, snpType = 'CasavaSNP')

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
					chromosome.set(number = currChrNumber, genome = genome, dataType = 'CasavaSNP')

				snp = CasavaSNP()
				snp.chromosome = chromosome
				snp.genome = genome
				#first column: chro, second first of range (identical to second column)
				snp.pos = int(sl[2])
				snp.bcalls_used = sl[3]
				snp.bcalls_filt = sl[4]
				snp.ref = sl[5]
				snp.QSNP = int(sl[6])
				snp.max_gt = uf.getPolymorphicNucleotide(sl[7])
				snp.Qmax_gt = int(sl[8])
				snp.max_gt_poly_site = sl[9]
				snp.Qmax_gt_poly_site = int(sl[10])
				snp.A_used = int(sl[11])
				snp.C_used = int(sl[12])
				snp.G_used = int(sl[13])
				snp.T_used = int(sl[14])
				snp.save()
				#genome.snps.append(snp)

		print 'saving...'
		genome.save()
		rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).commit()
		rabaDB.setup.RabaConnection(conf.pyGeno_RABA_NAMESPACE).setAutoCommit(True)
		print genome.snps
		print 'importation %s of genome %s done.' %(specie, genomeName)

	except (KeyboardInterrupt, SystemExit, Exception) :
		print "===>Exception caught! Rollback!<==="
		os.remove(conf.pyGeno_RABA_DBFILE)
		os.rename(bckFn, conf.pyGeno_RABA_DBFILE)
		raise

if __name__ == "__main__" :
	#deleteBackUps(forceDelete = True)(, )
	#deleteGenome(specie = 'Mus_musculus', name = 'GRCm38_test')
	#importGenome('/u/daoudat/py/pyGeno/importationPackages/genomes/mouse/mus_musculus_Y-only.tar.gz', verbose = 1)
	importGenome('/u/daoudat/py/pyGeno/importationPackages/GRCh37.74/GRCh37.74.tar.gz', verbose = 0)
	#g = Genome(specie = 'Mus_musculus', name = 'GRCm38_test')
	#print g.get(Transcript)[0].exons
	#deleteGenome(specie = 'human', name = 'ARN_R')
	#importGenome_casava('human', 'ARN_R', '/u/corona/Project_DSP008a/Build_Diana_ARN_R/snps.txt')
	#g = Genome(specie = 'human', name = 'ARN_R')
	#print g.snps
	#print Genome(specie = 'human', name = 'GRCh37.74')

