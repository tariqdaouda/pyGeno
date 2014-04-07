import os, glob, gzip, tarfile, shutil, time, sys, gc, cPickle, tempfile
from ConfigParser import SafeConfigParser

import configuration as conf
"""
import rabaDB.setup
from rabaDB.Raba import *
import rabaDB.fields as rf
from rabaDB.filters import RabaQuery
"""
from Genome import *
from Chromosome import *
from Gene import *
from Transcript import *
from Exon import *
from Protein import *

from pyGeno.tools.GTFTools import GTFFile
from pyGeno.tools.ProgressBar import ProgressBar
from pyGeno.tools.io import printf

def backUpDB() :
	"backup the current database version. automatically called by importGenome(). Returns the filename of the backup"
	st = time.ctime().replace(' ', '_')
	fn = conf.pyGeno_RABA_DBFILE.replace('.db', '_%s_auto-bck.db' % st)
	shutil.copy2(conf.pyGeno_RABA_DBFILE, fn)

	return fn


def _decompressPackage(packageFile) :
	pFile = tarfile.open(packageFile)
	packageDir = os.path.normpath('%s/tmp_pygeno_import' % tempfile.gettempdir())

 	if os.path.isdir(packageDir) :
		shutil.rmtree(packageDir)
	os.makedirs(packageDir)

	for mem in pFile :
		pFile.extract(mem, packageDir)

	return packageDir

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
		raise ValueError("There seems to be already a genome (%s, %s), please call deleteGenome() first if you want to reinstall it" % (genomeName, specie))
	except KeyError:
		pass

	seqTargetDir = conf.getGenomeSequencePath(specie.lower(), genomeName)
	if os.path.isdir(seqTargetDir) :
		raise ValueError("The directory %s already exists, Please call deleteGenome() first if you want to reinstall" % seqTargetDir)

	os.makedirs(seqTargetDir)

	genome = Genome_Raba()
	genome.set(name = genomeName, specie = specie, source = genomeSource, packageInfos = packageInfos)

	printf("Importing:\n\t%s\nGenome:\n\t%s\n..."  % (reformatItems(packageInfos).replace('\n', '\n\t'), reformatItems(parser.items('genome')).replace('\n', '\n\t')))

	chros = _importGenomeObjects(os.path.normpath(packageDir+'/'+gtfFile), chromosomeSet, genome, verbose)
	startChro = 0
	for chro in chros :
		printf("Importing DNA sequence of chromosome %s..." % chro)
		length = _importSequence(chro, os.path.normpath(packageDir+'/'+chromosomesFiles[chro.number.lower()]), seqTargetDir)
		chro.start = startChro
		chro.end = startChro+length
		startChro = chro.end

	genome.save()
	shutil.rmtree(packageDir)

def deleteGenome(specie, name) :
	"removes all infos about a genome"

	printf('deleting genome (%s, %s)...' % (specie, name))

	conf.db.beginTransaction()
	printf('\tdeleting genome information (%s, %s)...' % (specie, name))
	objs = []
	try :
		genome = Genome_Raba(name = name, specie = specie.lower())
		objs.append(genome)
		for typ in (Chromosome_Raba, Gene_Raba, Transcript_Raba, Exon_Raba, Protein_Raba) :
			f = RabaQuery(typ, namespace = genome._raba_namespace)
			f.addFilter({'genome' : genome})
			for e in f.iterRun() :
				objs.append(e)

		printf('\tdeleting...')
		for e in objs :
			e.delete()
	except KeyError as e :
		printf("\tWARNING, couldn't remove genome form db, maybe it's not there: ", e)

	printf('\tdeleting folder')
	try :
		shutil.rmtree(conf.getGenomeSequencePath(specie, name))
	except OSError as e:
		printf('\tWARNING, Unable to delete folder: ', e)

	conf.db.endTransaction()

def _importGenomeObjects(gtfFilePath, chroSet, genome, verbose = 0) :
	"verbose is int [0, 4] for various levels of verbosity"
	
	printf('Importing gene set infos from %s...' % gtfFilePath)
	startTime = time.time()
	
	printf('Backuping indexes...')
	indexes = conf.db.getIndexes()
	printf('Removing all indexes...')
	conf.db.flushIndexes()
	
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
				chromosomes[chroNumber] = Chromosome_Raba()
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
			frame = gtf.get(i, 'frame')

			start = int(gtf.get(i, 'start')) - 1
			end = int(gtf.get(i, 'end'))
			if start > end :
				start, end = end, start

			if geneId not in genes :
				if verbose > 0 :
					printf('\tGene %s, %s...' % (geneId, geneName))
				genes[geneId] = Gene_Raba()
				genes[geneId].set(genome = genome, id = geneId, chromosome = chromosomes[chroNumber], name = geneName, strand = strand, biotype = gene_biotype)

			try :
				transId = gtf.get(i, 'transcript_id')
				transName = gtf.get(i, 'transcript_name')
			except KeyError :
				if verbose > 2 :
					printf('\t\tWarning: no transcript_id, name found in line %d' % i)

			if transId not in transcripts :
				if verbose > 1 :
					printf('\t\tTranscript %s, %s...' % (transId, transName))
				transcripts[transId] = Transcript_Raba()
				transcripts[transId].set(genome = genome, id = transId, chromosome = chromosomes[chroNumber], gene = genes[geneId], name = transName)
			try :
				protId = gtf.get(i, 'protein_id')
				if protId not in proteins :
					if verbose > 1 :
						printf('\t\tProtein %s...' % (protId))
					proteins[protId] = Protein_Raba()
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
						exons[exonKey] = Exon_Raba()
						exons[exonKey].set(genome = genome, id = exonId, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], protein = proteins[protId], strand = strand, number = exonNumber, start = start, end = end)
						transcripts[transId].exons.append(exons[exonKey])

				except KeyError :
					printf('Warning: no exon_id found in line %d' % i)

			elif regionType == 'CDS' :
				exons[exonKey].CDS_start = start
				exons[exonKey].CDS_end = end
				exons[exonKey].frame = frame
			elif regionType == 'stop_codon' :
				if strand == '+' :
					exons[exonKey].end += 3
					if exons[exonKey].CDS_end != None :
						exons[exonKey].CDS_end += 3
				if strand == '-' :
					exons[exonKey].start -= 3
					if exons[exonKey].CDS_start != None :
						exons[exonKey].CDS_start -= 3
	
	printf('\tdone (%fmin), total time (%f).' %((time.time()-chroStartTime)/60, (time.time()-startTime)/60))

	conf.db.beginTransaction()
	printf('saving %d chromsomes...' % len(chromosomes))
	for c in chromosomes.itervalues() :
		if verbose > 1 :
			printf('saving chromsomes %s...' % c.number)
		c.save()

	printf('saving %d genes...' % len(genes))
	i = 0
	for c in genes.itervalues() :
		if verbose > 1 :
			printf('saving gene %s (%f%%)...' % (c.id, float(i)/len(transcripts)*100))
		i += 1
		c.save()

	printf('saving %d transcripts and exons %d...' % (len(transcripts), len(exons)))
	i = 0
	for c in transcripts.itervalues() :
		if verbose > 1 :
			printf('saving transcripts %s (%f%%)...' % (c.id, float(i)/len(transcripts)*100))
		i += 1
		c.save()
	conf.db.enableQueryPrint(False)

	printf('saving %d proteins...' % len(proteins))
	i = 0
	for c in proteins.itervalues() :
		if verbose > 1 :
			printf('saving protein %s (%f%%)...' % (c.id, float(i)/len(transcripts)*100))
		i += 1
		c.save()

	printf('saving genome...')
	genome.save()
	
	printf('restoring core indexes...')
	Transcript.ensureGlobalIndex('exons')
	Chromosome.ensureGlobalIndex('genome')
	Gene.ensureGlobalIndex('genome')
	Transcript.ensureGlobalIndex('genome')
	Protein.ensureGlobalIndex('genome')
	Exon.ensureGlobalIndex('genome')
	
	printf('restoring prior indexes')
	for idx in indexes :
		conf.db.execute(idx[-1])

	printf('\tcommiting changes...')
	conf.db.endTransaction()
	printf('\tdone total time (%f).' %((time.time()-startTime)/60))

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



if __name__ == "__main__" :
	print """
	Ex :
	deleteGenome(specie = 'Mus_musculus', name = 'GRCm38_test')
	importGenome('/u/daoudat/py/pyGeno/importationPackages/genomes/mouse/mus_musculus_Y-only.tar.gz', verbose = 0)"""
