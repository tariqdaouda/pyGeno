import os, glob, gzip, shutil, time
from ConfigParser import SafeConfigParser

import configuration as conf

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
	st = time.ctime().replace(' ', '_')
	fn = conf.pyGeno_RABA_DBFILE.replace('.db', '%s_auto-bck.db' % st)
	shutil.copy2(conf.pyGeno_RABA_DBFILE, fn)
	
	return fn

def importGenome(packageDir, verbose = False) :
	"""TODO: Write description!"""
	
	def reformatItems(items) :
		s = str(items)
		s = s.replace('[', '').replace(']', '').replace("',", ': ').replace('), ', '\n').replace("'", '').replace('(', '').replace(')', '')
		return s
	

	parser = SafeConfigParser()
	parser.read(packageDir+'/manifest.ini')
	packageInfos = parser.items('package_infos')
	
	genomeName = parser.get('genome', 'name')
	specie = parser.get('genome', 'specie')	
	genomeSource = parser.get('genome', 'source')
	gtfFile = parser.get('gene_set', 'gtf')
	chromosomesFiles = dict(parser.items('chromosome_files'))
	chromosomeSet = set(chromosomesFiles.keys())
	
	print "Installing:\n\t%s\nGenome:\n\t%s\n..."  % (reformatItems(packageInfos).replace('\n', '\n\t'), reformatItems(parser.items('genome')).replace('\n', '\n\t'))
	bckFn = backUpDB()
	print "=====\nIf anything goes wrong, the db has been backuped here: %s\nSimply rename it to: %s\n=====" %(bckFn, conf.pyGeno_RABA_DBFILE)

	genome = Genome()
	genome.set(name = genomeName, specie = specie, source = genomeSource, packageInfos = packageInfos)
	seqTargetDir = genome.getSequencePath()
	
	#if os.path.isdir(seqTargetDir) :
	#	raise ValueError("The directory %s already exists. If you want to reinstall a package delete the folder first" % seqTargetDir)
	#os.makedirs(seqTargetDir)
	
	_importGenomeObjects(packageDir+'/'+gtfFile, chromosomeSet, genome, verbose)
	x1Chro = 0
	for chro in genome.chromosomes :
		length = __importSequence(chro, packageDir+'/'+chromosomesFiles[chro.number.lower()], seqTargetDir)
		chro.x1 = x1Chro
		chro.x2 = x1Chro+length
		x1Chro = chro.x2
	
	genome.save()
	
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
		chroNumber = gtf.get(i, 'seqname')
		if chroNumber.upper() in chroSet or chroNumber.lower() in chroSet:
			
			chroNumber = chroNumber#.lower()
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
			x1 = int(gtf.get(i, 'start'))
			x2 = int(gtf.get(i, 'end')) +1

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
				transcripts[transId] = Transcript()
				transcripts[transId].set(genome = genome, id = transId, chromosome = chromosomes[chroNumber], gene = genes[geneId], name = transName)
			try :
				protId = gtf.get(i, 'protein_id')
				if protId not in proteins :
					if verbose :
						print '\tProtein %s...' % (protId)
					proteins[protId] = Protein()
					proteins[protId].set(genome = genome, id = protId, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], name = transName)
					transcripts[transId].protein = proteins[protId]
					proteins[protId].save()
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
						exons[exonKey] = Exon()
						exons[exonKey].set(genome = genome, id = exonId, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], strand = strand, number = exonNumber, x1 = x1, x2 = x2)
						exons[exonKey].save()
				except KeyError :
					print 'Warning: no exon_id found in line %d' % i
	
			elif regionType == 'CDS' :
				exons[exonKey].setCDS(x1, x2)
			elif regionType == 'start_codon' :
				exons[exonKey].startCodon = x1
			elif regionType == 'stop_codon' :
				exons[exonKey].stopCodon = x2
	
	print 'creating relations...'
	print '\ttranscript.exons...'
	for transcript in transcripts.values() :
		f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, Exon)
		f.addFilter(**{'transcript' : transcript})
		transcript.exons = f.run()
	
	print '\tgene.transcripts...'
	for gene in genes.values() :
		f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, Transcript)
		f.addFilter(**{'gene' : gene})
		gene.transcripts = f.run()
		
		f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, Exon)
		f.addFilter(**{'gene' : gene})
		gene.exons = f.run()
	
	print '\tchromosome.genes...'
	for chro in chromosomes.values() :
		f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, Gene)
		f.addFilter(**{'chromosome' : chro})
		chro.genes = f.run()
	
	print '\tgenome.chromosomes...'
	genome.chromosomes = RabaList(chromosomes.values())
	print 'saving...'
	print 'Done.'
	genome.save()
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
				chromosome.dataType = 'light'
				
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
		- all SNPs that have a chromosome different from the one of the file, or chromosome number of '?' are discarded
		- the default value for numeric values (ex: het, se(het), 'MAF', 'GMAf count') is 0.0
		- no '?' allowed.If a numeric has a value of '?' (ex: het, se(het), 'MAF', 'GMAf count') the value is set to 0.0
		- all snps have an orientation of +. If a snp has an orientation of -, it's alleles are replaced by their complements
		- positions are 0 based
		- the only value extracted from the files are : 'posistion', 'rs', 'type', 'assembly', 'chromosome', 'validated', 'alleles', 'original_orientation', 'maf_allele', 'maf_count', 'maf', 'het', 'se(het)
		-loc is a dictionary allele wise that simplifies the line loc' 
	"""
	
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
					pass
				
			elif sl[0][:4] == 'GMAF' :
				snp.maf_allele = sl[1].strip().replace('allele=', '')
				maf_count = sl[2].strip().replace('count=', '')
				try :
					snp.maf_count = float(maf_count)
				except :
					pass
				
				maf = sl[3].strip().replace('MAF=', '')
				try :
					snp.maf = float(maf)
				except :
					pass
			
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

	files = glob.glob(packageFolder+'/*.flat.gz')

	for fil in files :
		chrStrStartPos = fil.find('ch')
		chrStrStopPos = fil.find('.flat')
		numericFieldsWithNonNumericValues = 0
		
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
	importGenome('/u/daoudat/py/pyGeno/importationPackages/genomes/mouse/mus_muslcus/', verbose = False)
	g = Genome(specie = 'Mus_musculus', name = 'GRCm38_test')
	a = g.get(Exon, {'x1 >' : 797276})
	#print a[0]
	#importGenome_casava('Mus_musculus', 'musCasava', 'http://www.bioinfo.iric.ca/seq/Project_DSP008a/Build_Diana_ARN_R/snps.txt')
	
	#import_dbSNP('/u/daoudat/py/pyGeno/importationPackages/dbSNP/human/test', 'Mus_musculus', 'test')
	#print g.chromosomes
	#f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, dbSNP_SNP)
	#f.addFilter(**{'version' : 'test'})
	#for e in f.run() :
	#	print e.rsId
