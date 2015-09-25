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

import gc
#~ import objgraph

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

def deleteGenome(species, name) :
	"""Removes a genome from the database"""

	printf('deleting genome (%s, %s)...' % (species, name))

	conf.db.beginTransaction()
	objs = []
	allGood = True
	try :
		genome = Genome_Raba(name = name, species = species.lower())
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
		shutil.rmtree(conf.getGenomeSequencePath(species, name))
	except OSError as e:
		#~ printf('\tWARNING, Unable to delete folder: ', e)
		OSError('\tWARNING, Unable to delete folder: ', e)
		allGood = False
		
	conf.db.endTransaction()
	return allGood

def importGenome(packageFile, batchSize = 50, verbose = 0) :
	"""Import a pyGeno genome package. A genome packages is a tar.gz ball that contains at it's root:

	* gziped fasta files for all chromosomes, or URLs from where them must be downloaded
	
	* gziped GTF gene_set file from Ensembl, or an URL from where it must be downloaded
	
	* a manifest.ini file such as::
	
		[package_infos]
		description = Test package. This package installs only chromosome Y of mus musculus
		maintainer = Tariq Daouda
		maintainer_contact = tariq.daouda [at] umontreal
		version = GRCm38.73

		[genome]
		species = Mus_musculus
		name = GRCm38_test
		source = http://useast.ensembl.org/info/data/ftp/index.html

		[chromosome_files]
		Y = Mus_musculus.GRCm38.73.dna.chromosome.Y.fa.gz / or an url such as ftp://... or http://

		[gene_set]
		gtf = Mus_musculus.GRCm38.73_Y-only.gtf.gz / or an url such as ftp://... or http://

	All files except the manifest can be downloaded from: http://useast.ensembl.org/info/data/ftp/index.html
	
	A rollback is performed if an exception is caught during importation
	
	batchSize sets the number of genes to parse before performing a database save. PCs with little ram like
	small values, while those endowed with more memory may perform faster with higher ones.
	
	Verbose must be an int [0, 4] for various levels of verbosity
	"""

	def reformatItems(items) :
		s = str(items)
		s = s.replace('[', '').replace(']', '').replace("',", ': ').replace('), ', '\n').replace("'", '').replace('(', '').replace(')', '')
		return s

	printf('Importing genome package: %s... (This may take a while)' % packageFile)

	packageDir = _decompressPackage(packageFile)

	parser = SafeConfigParser()
	parser.read(os.path.normpath(packageDir+'/manifest.ini'))
	packageInfos = parser.items('package_infos')

	genomeName = parser.get('genome', 'name')
	species = parser.get('genome', 'species')
	genomeSource = parser.get('genome', 'source')
	
	seqTargetDir = conf.getGenomeSequencePath(species.lower(), genomeName)
	if os.path.isdir(seqTargetDir) :
		raise KeyError("The directory %s already exists, Please call deleteGenome() first if you want to reinstall" % seqTargetDir)
		
	gtfFile = _getFile(parser.get('gene_set', 'gtf'), packageDir)
	
	chromosomesFiles = {}
	chromosomeSet = set()
	for key, fil in parser.items('chromosome_files') :
		chromosomesFiles[key] = _getFile(fil, packageDir)
		chromosomeSet.add(key)

	try :
		genome = Genome(name = genomeName, species = species)
	except KeyError:
		pass
	else :
		raise KeyError("There seems to be already a genome (%s, %s), please call deleteGenome() first if you want to reinstall it" % (genomeName, species))

	genome = Genome_Raba()
	genome.set(name = genomeName, species = species, source = genomeSource, packageInfos = packageInfos)

	printf("Importing:\n\t%s\nGenome:\n\t%s\n..."  % (reformatItems(packageInfos).replace('\n', '\n\t'), reformatItems(parser.items('genome')).replace('\n', '\n\t')))

	chros = _importGenomeObjects(gtfFile, chromosomeSet, genome, batchSize, verbose)
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
	
	#~ objgraph.show_most_common_types(limit=20)
	return True

#~ @profile
def _importGenomeObjects(gtfFilePath, chroSet, genome, batchSize, verbose = 0) :
	"""verbose must be an int [0, 4] for various levels of verbosity"""

	class Store(object) :
		
		def __init__(self, conf) :
			self.conf = conf
			self.chromosomes = {}
			
			self.genes = {}
			self.transcripts = {}
			self.proteins = {}
			self.exons = {}

		def batch_save(self) :
			self.conf.db.beginTransaction()
			
			for c in self.genes.itervalues() :
				c.save()
				conf.removeFromDBRegistery(c)
				
			for c in self.transcripts.itervalues() :
				c.save()
				conf.removeFromDBRegistery(c.exons)
				conf.removeFromDBRegistery(c)
			
			for c in self.proteins.itervalues() :
				c.save()
				conf.removeFromDBRegistery(c)
			
			self.conf.db.endTransaction()
			
			del(self.genes)
			del(self.transcripts)
			del(self.proteins)
			del(self.exons)
			
			self.genes = {}
			self.transcripts = {}
			self.proteins = {}
			self.exons = {}

			gc.collect()

		def save_chros(self) :
			pBar = ProgressBar(nbEpochs = len(self.chromosomes))
			for c in self.chromosomes.itervalues() :
				pBar.update(label = 'Chr %s' % c.number)
				c.save()
			pBar.close()
		
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
	
	store = Store(conf)
	chroNumber = None
	pBar = ProgressBar(nbEpochs = len(gtf))
	for line in gtf :
		chroN = line['seqname']
		pBar.update(label = "Chr %s" % chroN)
		
		if (chroN.upper() in chroSet or chroN.lower() in chroSet):
			strand = line['strand']
			gene_biotype = line['gene_biotype']
			regionType = line['feature']
			frame = line['frame']

			start = int(line['start']) - 1
			end = int(line['end'])
			if start > end :
				start, end = end, start

			chroNumber = chroN.upper()
			if chroNumber not in store.chromosomes :
				store.chromosomes[chroNumber] = Chromosome_Raba()
				store.chromosomes[chroNumber].set(genome = genome, number = chroNumber)
			
			try :
				geneId = line['gene_id']
				geneName =  line['gene_name']
			except KeyError :
				geneId = None
				geneName = None
				if verbose :
					printf('Warning: no gene_id/name found in line %s' % gtf[i])
			
			if geneId is not None :
				if geneId not in store.genes :
					if len(store.genes) > batchSize :
						store.batch_save()
					
					if verbose > 0 :
						printf('\tGene %s, %s...' % (geneId, geneName))
					store.genes[geneId] = Gene_Raba()
					store.genes[geneId].set(genome = genome, id = geneId, chromosome = store.chromosomes[chroNumber], name = geneName, strand = strand, biotype = gene_biotype)
				if start < store.genes[geneId].start or store.genes[geneId].start is None :
						store.genes[geneId].start = start
				if end > store.genes[geneId].end or store.genes[geneId].end is None :
					store.genes[geneId].end = end
			try :
				transId = line['transcript_id']
				transName = line['transcript_name']
			except KeyError :
				transId = None
				transName = None
				if verbose > 2 :
					printf('\t\tWarning: no transcript_id, name found in line %s' % gtf[i])
			
			if transId is not None :
				if transId not in store.transcripts :
					if verbose > 1 :
						printf('\t\tTranscript %s, %s...' % (transId, transName))
					store.transcripts[transId] = Transcript_Raba()
					store.transcripts[transId].set(genome = genome, id = transId, chromosome = store.chromosomes[chroNumber], gene = store.genes.get(geneId, None), name = transName)
				if start < store.transcripts[transId].start or store.transcripts[transId].start is None:
					store.transcripts[transId].start = start
				if end > store.transcripts[transId].end or store.transcripts[transId].end is None:
					store.transcripts[transId].end = end
			
				try :
					protId = line['protein_id']
				except KeyError :
					protId = None
					if verbose > 2 :
						printf('Warning: no protein_id found in line %s' % gtf[i])
				
				if protId is not None and protId not in store.proteins :
					if verbose > 1 :
						printf('\t\tProtein %s...' % (protId))
					store.proteins[protId] = Protein_Raba()
					store.proteins[protId].set(genome = genome, id = protId, chromosome = store.chromosomes[chroNumber], gene = store.genes.get(geneId, None), transcript = store.transcripts.get(transId, None), name = transName)
					store.transcripts[transId].protein = store.proteins[protId]

				try :
					exonNumber = int(line['exon_number']) - 1
					exonKey = (transId, exonNumber)
				except KeyError :
					exonNumber = None
					exonKey = None
					if verbose > 2 :
						printf('Warning: no exon number or id found in line %s' % gtf[i])
				
				if exonKey is not None :
					if verbose > 3 :
						printf('\t\t\texon %s...' % (exonId))
					
					if exonKey not in store.exons :
						store.exons[exonKey] = Exon_Raba()
						store.exons[exonKey].set(genome = genome, chromosome = store.chromosomes[chroNumber], gene = store.genes.get(geneId, None), transcript = store.transcripts.get(transId, None), protein = store.proteins.get(protId, None), strand = strand, number = exonNumber, start = start, end = end)
						store.transcripts[transId].exons.append(store.exons[exonKey])
					
					try :
						store.exons[exonKey].id = line['exon_id']
					except KeyError :
						pass
					
					if regionType == 'exon' :
						if start < store.exons[exonKey].start or store.exons[exonKey].start is None:
							store.exons[exonKey].start = start
						if end > store.transcripts[transId].end or store.exons[exonKey].end is None:
							store.exons[exonKey].end = end
					elif regionType == 'CDS' :
						store.exons[exonKey].CDS_start = start
						store.exons[exonKey].CDS_end = end
						store.exons[exonKey].frame = frame
					elif regionType == 'stop_codon' :
						if strand == '+' :
							store.exons[exonKey].end += 3
							if store.exons[exonKey].CDS_end != None :
								store.exons[exonKey].CDS_end += 3
						if strand == '-' :
							if store.exons[exonKey].CDS_start != None :
								store.exons[exonKey].CDS_start -= 3
	pBar.close()
	
	store.batch_save()
	
	conf.db.beginTransaction()
	printf('almost done saving chromosomes...')
	store.save_chros()
	
	printf('saving genome object...')
	genome.save()
	conf.db.endTransaction()
	
	conf.db.beginTransaction()
	printf('restoring core indexes...')
	# Genome.ensureGlobalIndex(('name', 'species'))
	# Chromosome.ensureGlobalIndex('genome')
	# Gene.ensureGlobalIndex('genome')
	# Transcript.ensureGlobalIndex('genome')
	# Protein.ensureGlobalIndex('genome')
	# Exon.ensureGlobalIndex('genome')
	Transcript.ensureGlobalIndex('exons')
	
	printf('commiting changes...')
	conf.db.endTransaction()
	
	conf.db.beginTransaction()
	printf('restoring user indexes')
	pBar = ProgressBar(label = "restoring", nbEpochs = len(indexes))
	for idx in indexes :
		pBar.update()
		conf.db.execute(idx[-1].replace('CREATE INDEX', 'CREATE INDEX IF NOT EXISTS'))
	pBar.close()
	
	printf('commiting changes...')
	conf.db.endTransaction()
	
	return store.chromosomes.values()

#~ @profile
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
