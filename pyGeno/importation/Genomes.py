import os, glob, gzip, tarfile, shutil, time, sys, gc, cPickle, tempfile, urllib
from ConfigParser import SafeConfigParser

from pyGeno.tools.ProgressBar import ProgressBar
import pyGeno.configuration as conf

from pyGeno.Genome import *
from pyGeno.Chromosome import *
from pyGeno.Gene import *
from pyGeno.Transcript import *
from pyGeno.Exon import *
from pyGeno.Protein import *

from pyGeno.tools.parsers.GTFTools import GTFFile
from pyGeno.tools.ProgressBar import ProgressBar
from pyGeno.tools.io import printf

def backUpDB() :
	"""backup the current database version. automatically called by importGenome(). Returns the filename of the backup"""
	st = time.ctime().replace(' ', '_')
	fn = conf.pyGeno_RABA_DBFILE.replace('.db', '_%s-bck.db' % st)
	shutil.copy2(conf.pyGeno_RABA_DBFILE, fn)

	return fn

def _decompressPackage(packageFile) :
	pFile = tarfile.open(packageFile)
	
	packageDir = tempfile.mkdtemp(prefix = "pyGeno_import_")
 	if os.path.isdir(packageDir) :
		shutil.rmtree(packageDir)
	os.makedirs(packageDir)

	for mem in pFile :
		pFile.extract(mem, packageDir)

	return packageDir

def _getFile(fil, directory) :
	if fil.find("http://") == 0 or fil.find("ftp://") == 0 :
		printf("Downloading file: %s..." % fil)
		finalFile = os.path.normpath('%s/%s' %(directory, fil.split('/')[-1]))
		urllib.urlretrieve (fil, finalFile)
		printf('done.')
	else :
		finalFile = os.path.normpath('%s/%s' %(directory, fil))
	
	return finalFile

def deleteGenome(specie, name) :
	"""Removes a genome from the database"""

	printf('deleting genome (%s, %s)...' % (specie, name))

	conf.db.beginTransaction()
	objs = []
	allGood = True
	try :
		genome = Genome_Raba(name = name, specie = specie.lower())
		objs.append(genome)
		pBar = ProgressBar(label = 'preparing')
		for typ in (Chromosome_Raba, Gene_Raba, Transcript_Raba, Exon_Raba, Protein_Raba) :
			pBar.update()
			f = RabaQuery(typ, namespace = genome._raba_namespace)
			f.addFilter({'genome' : genome})
			for e in f.iterRun() :
				objs.append(e)
		pBar.close()
		
		pBar = ProgressBar(nbEpochs = len(objs), label = 'deleting objects')
		for e in objs :
			pBar.update()
			e.delete()
		pBar.close()
		
	except KeyError as e :
		#~ printf("\tWARNING, couldn't remove genome form db, maybe it's not there: ", e)
		raise KeyError("\tWARNING, couldn't remove genome form db, maybe it's not there: ", e)
		allGood = False
	printf('\tdeleting folder')
	try :
		shutil.rmtree(conf.getGenomeSequencePath(specie, name))
	except OSError as e:
		#~ printf('\tWARNING, Unable to delete folder: ', e)
		OSError('\tWARNING, Unable to delete folder: ', e)
		allGood = False
		
	conf.db.endTransaction()
	return allGood

def importGenome(packageFile, verbose = False) :
	"""Import a pyGeno genome package. A genome packages is a tar.gz ball that contains at it's root:

	* gziped fasta files for all chromosomes
	
	* gziped GTF gene_set file from ensembl
	
	* a manifest.ini file such as::
	
		[package_infos]
		description = Test package. This package installs only chromosome Y of mus musculus
		maintainer = Tariq Daouda
		maintainer_contact = tariq.daouda [at] umontreal
		version = GRCm38.73

		[genome]
		specie = Mus_musculus
		name = GRCm38_test
		source = http://useast.ensembl.org/info/data/ftp/index.html

		[chromosome_files]
		Y = Mus_musculus.GRCm38.73.dna.chromosome.Y.fa.gz / or an url such as ftp://... or http://

		[gene_set]
		gtf = Mus_musculus.GRCm38.73_Y-only.gtf.gz / or an url such as ftp://... or http://

	All files except the manifest can be downloaded from: http://useast.ensembl.org/info/data/ftp/index.html
	
	A rollback is performed if an exception is caught during importation
	"""

	def reformatItems(items) :
		s = str(items)
		s = s.replace('[', '').replace(']', '').replace("',", ': ').replace('), ', '\n').replace("'", '').replace('(', '').replace(')', '')
		return s

	printf('Importing genome package: %s... (This make take a while)' % packageFile)

	packageDir = _decompressPackage(packageFile)

	parser = SafeConfigParser()
	parser.read(os.path.normpath(packageDir+'/manifest.ini'))
	packageInfos = parser.items('package_infos')

	genomeName = parser.get('genome', 'name')
	specie = parser.get('genome', 'specie')
	genomeSource = parser.get('genome', 'source')
	
	seqTargetDir = conf.getGenomeSequencePath(specie.lower(), genomeName)
	if os.path.isdir(seqTargetDir) :
		raise ValueError("The directory %s already exists, Please call deleteGenome() first if you want to reinstall" % seqTargetDir)
		
	gtfFile = _getFile(parser.get('gene_set', 'gtf'), packageDir)
	
	chromosomesFiles = {}
	chromosomeSet = set()
	for key, fil in parser.items('chromosome_files') :
		chromosomesFiles[key] = _getFile(fil, packageDir)
		chromosomeSet.add(key)

	try :
		genome = Genome(name = genomeName, specie = specie)
		raise ValueError("There seems to be already a genome (%s, %s), please call deleteGenome() first if you want to reinstall it" % (genomeName, specie))
	except KeyError:
		pass

	genome = Genome_Raba()
	genome.set(name = genomeName, specie = specie, source = genomeSource, packageInfos = packageInfos)

	printf("Importing:\n\t%s\nGenome:\n\t%s\n..."  % (reformatItems(packageInfos).replace('\n', '\n\t'), reformatItems(parser.items('genome')).replace('\n', '\n\t')))

	chros = _importGenomeObjects(gtfFile, chromosomeSet, genome, verbose)
	os.makedirs(seqTargetDir)
	startChro = 0
	pBar = ProgressBar(nbEpochs = len(chros))
	for chro in chros :
		pBar.update(label = "Importing DNA, chro %s" % chro.number)
		length = _importSequence(chro, chromosomesFiles[chro.number.lower()], seqTargetDir)
		chro.start = startChro
		chro.end = startChro+length
		startChro = chro.end
	pBar.close()
	
	shutil.rmtree(packageDir)
	return True
	
def _importGenomeObjects(gtfFilePath, chroSet, genome, verbose = 0) :
	"verbose is int [0, 4] for various levels of verbosity"
	
	printf('Importing gene set infos from %s...' % gtfFilePath)
	
	printf('Backuping indexes...')
	indexes = conf.db.getIndexes()
	printf("Droping all your indexes, (don't worry i'll restore them later)...")
	Genome_Raba.flushIndexes()
	Chromosome_Raba.flushIndexes()
	Gene_Raba.flushIndexes()
	Transcript_Raba.flushIndexes()
	Protein_Raba.flushIndexes()
	Exon_Raba.flushIndexes()
	
	printf("Parsing gene set...")
	gtf = GTFFile(gtfFilePath, gziped = True)
	printf('Done. Importation begins!')
	
	chromosomes = {}
	genes = {}
	transcripts = {}
	proteins = {}
	exons = {}
	chroNumber = None
	
	pBar = ProgressBar(nbEpochs = len(gtf))
	for i in xrange(len(gtf)) :
		chroN = str(gtf.get(i, 'seqname'))
		pBar.update(label = "Chr %s" % chroN)
		
		if (chroN.upper() in chroSet or chroN.lower() in chroSet):
			strand = gtf.get(i, 'strand')
			gene_biotype = gtf.get(i, 'gene_biotype')
			regionType = gtf.get(i, 'feature')
			frame = gtf.get(i, 'frame')

			start = int(gtf.get(i, 'start')) - 1
			end = int(gtf.get(i, 'end'))
			if start > end :
				start, end = end, start

			chroNumber = chroN.upper()
			if chroNumber not in chromosomes :
				chromosomes[chroNumber] = Chromosome_Raba()
				chromosomes[chroNumber].set(genome = genome, number = chroNumber)
				chromosomes[chroNumber].dataType = 'heavy'

			try :
				geneId = gtf.get(i, 'gene_id')
				geneName = gtf.get(i, 'gene_name')
			except KeyError :
				geneId = None
				geneName = None
				if verbose :
					printf('Warning: no gene_id/name found in line %s' % gtf[i])
			
			if geneId is not None :
				if geneId not in genes :
					if verbose > 0 :
						printf('\tGene %s, %s...' % (geneId, geneName))
					genes[geneId] = Gene_Raba()
					genes[geneId].set(genome = genome, id = geneId, chromosome = chromosomes[chroNumber], name = geneName, strand = strand, biotype = gene_biotype)
				if start < genes[geneId].start or genes[geneId].start is None :
						genes[geneId].start = start
				if end > genes[geneId].end or genes[geneId].end is None :
					genes[geneId].end = end
			try :
				transId = gtf.get(i, 'transcript_id')
				transName = gtf.get(i, 'transcript_name')
			except KeyError :
				transId = None
				transName = None
				if verbose > 2 :
					printf('\t\tWarning: no transcript_id, name found in line %s' % gtf[i])
			
			if transId is not None :
				if transId not in transcripts :
					if verbose > 1 :
						printf('\t\tTranscript %s, %s...' % (transId, transName))
					transcripts[transId] = Transcript_Raba()
					transcripts[transId].set(genome = genome, id = transId, chromosome = chromosomes[chroNumber], gene = genes.get(geneId, None), name = transName)
				if start < transcripts[transId].start or transcripts[transId].start is None:
					transcripts[transId].start = start
				if end > transcripts[transId].end or transcripts[transId].end is None:
					transcripts[transId].end = end
			
				try :
					protId = gtf.get(i, 'protein_id')
				except KeyError :
					protId = None
					if verbose > 2 :
						printf('Warning: no protein_id found in line %s' % gtf[i])
				
				if protId is not None and protId not in proteins :
					if verbose > 1 :
						printf('\t\tProtein %s...' % (protId))
					proteins[protId] = Protein_Raba()
					proteins[protId].set(genome = genome, id = protId, chromosome = chromosomes[chroNumber], gene = genes.get(geneId, None), transcript = transcripts.get(transId, None), name = transName)
					transcripts[transId].protein = proteins[protId]

				try :
					exonNumber = int(gtf.get(i, 'exon_number')) - 1
					exonKey = (transId, exonNumber)
				except KeyError :
					exonNumber = None
					exonKey = None
					if verbose > 2 :
						printf('Warning: no exon number or id found in line %s' % gtf[i])
				
				if exonKey is not None :
					if verbose > 3 :
						printf('\t\t\texon %s...' % (exonId))
					
					if exonKey not in exons :
						exons[exonKey] = Exon_Raba()
						exons[exonKey].set(genome = genome, chromosome = chromosomes[chroNumber], gene = genes.get(geneId, None), transcript = transcripts.get(transId, None), protein = proteins.get(protId, None), strand = strand, number = exonNumber, start = start, end = end)
						transcripts[transId].exons.append(exons[exonKey])
					
					try :
						exons[exonKey].id = gtf.get(i, 'exon_id')
					except KeyError :
						pass
					
					if regionType == 'exon' :
						if start < exons[exonKey].start or exons[exonKey].start is None:
							exons[exonKey].start = start
						if end > transcripts[transId].end or exons[exonKey].end is None:
							exons[exonKey].end = end
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
							if exons[exonKey].CDS_start != None :
								exons[exonKey].CDS_start -= 3
	
	pBar.close()
	conf.db.beginTransaction()
	printf('saving %d chromsomes...' % len(chromosomes))
	pBar = ProgressBar(nbEpochs = len(chromosomes))
	for c in chromosomes.itervalues() :
		pBar.update(label = 'Chr %s' % c.number)
		c.save()
	pBar.close()
	
	printf('saving %d genes...' % len(genes))
	pBar = ProgressBar(nbEpochs = len(genes))
	for c in genes.itervalues() :
		pBar.update(label = 'gene %s' % c.id)
		c.save()
	pBar.close()

	printf('saving %d transcripts and exons %d...' % (len(transcripts), len(exons)))
	pBar = ProgressBar(nbEpochs = len(transcripts))
	for c in transcripts.itervalues() :
		pBar.update(label = 'transcript %s' % c.id)
		c.save()
	pBar.close()
	
	printf('saving %d proteins...' % len(proteins))
	pBar = ProgressBar(nbEpochs = len(proteins))
	for c in proteins.itervalues() :
		pBar.update(label = 'protein %s' % c.id)
		c.save()
	pBar.close()

	printf('saving genome...')
	genome.save()
	
	printf('restoring core indexes...')
	Genome.ensureGlobalIndex(('name', 'specie'))
	Chromosome.ensureGlobalIndex('genome')
	Gene.ensureGlobalIndex('genome')
	Transcript.ensureGlobalIndex('genome')
	Transcript.ensureGlobalIndex('exons')
	Protein.ensureGlobalIndex('genome')
	Exon.ensureGlobalIndex('genome')
	
	printf('commiting changes...')
	conf.db.endTransaction()
	
	conf.db.beginTransaction()
	printf('restoring user indexes')
	pBar = ProgressBar(nbEpochs = len(indexes))
	for idx in indexes :
		pBar.update()
		conf.db.execute(idx[-1].replace('CREATE INDEX', 'CREATE INDEX IF NOT EXISTS'))
	pBar.close()
	
	printf('commiting changes...')
	conf.db.endTransaction()
	
	return chromosomes.values()

def _importSequence(chromosome, fastaFile, targetDir) :
	"Serializes fastas into .dat files"

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
