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

#from tools import UsefulFunctions as uf
#from tools.CSVTools import CSVFile

from tools.GTFTools import GTFFile

#from expyutils.GTFTools import GTFFile

def backUpDB() :
	st = time.ctime().replace(' ', '_')
	fn = conf.pyGeno_RABA_DBFILE.replace('.db', '%s_auto-bck.db' % st)
	shutil.copy2(conf.pyGeno_RABA_DBFILE, fn)
	
	return fn

#===These are the only function that you will need to import new features  in pyGeno===

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

	genome = Genome(name = genomeName, specie = specie, genomeSource = genomeSource)
	genome.packageInfos = packageInfos
	seqTargetDir = genome.getSequencePath()
	
	#if os.path.isdir(seqTargetDir) :
	#	raise ValueError("The directory %s already exists. If you want to reinstall a package delete the folder first" % seqTargetDir)
	
	#os.makedirs(seqTargetDir)
	
	_importGenomeObjects(packageDir+'/'+gtfFile, chromosomeSet, genome, verbose)
	#print genome.chromosomes
	x1Chro = 0
	for chro in genome.chromosomes :
		#print chro
		length = __importSequence(chro, packageDir+'/'+chromosomesFiles[chro.number.lower()], seqTargetDir)
		chro.x1 = x1Chro
		chro.x2 = x1Chro+length
		x1Chro = chro.x2
		#chro.save()
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
				chromosomes[chroNumber] = Chromosome(genome = genome, number = chroNumber)
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
				genes[geneId] = Gene(genome = genome, id = geneId)
				genes[geneId].set(chromosome = chromosomes[chroNumber], name = geneName, strand = strand, biotype = gene_biotype)
			
			try :
				transId = gtf.get(i, 'transcript_id')
				transName = gtf.get(i, 'transcript_name')
			except KeyError :
				if verbose :
					print '\tWarning: no transcript_id, name found in line %d' % i
			
			if transId not in transcripts :
				if verbose :
					print '\tTranscript %s, %s...' % (transId, transName)
				transcripts[transId] = Transcript(genome = genome, id = transId)
				##print transcripts[transId]
				transcripts[transId].set(chromosome = chromosomes[chroNumber], gene = genes[geneId], name = transName)
				#print transcripts[transId].exons
			try :
				protId = gtf.get(i, 'protein_id')
				if protId not in proteins :
					if verbose :
						print '\tProtein %s...' % (protId)
					proteins[protId] = Protein(genome = genome, id = protId)
					proteins[protId].set(chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], name = transName)
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
						exons[exonKey] = Exon(genome = genome, id = exonId)
						exons[exonKey].set(chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], strand = strand, number = exonNumber, x1 = x1, x2 = x2)
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
	genome.save()
	print 'Done.'
	
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
	"""Creates a light genome (contains only snps infos and no sequence from the reference genome)
	The .casavasnps files generated are identical to the casava snps but with ';' instead of tabs and 
	a single position instead of a range"""

	path = conf.pyGeno_SETTINGS['DATA_PATH']+'/%s/genomes/%s/'%(specie, genomeName)
	print 'importing genome %s...' %path
	
	if not os.path.exists(path):
		os.makedirs(path)
	
	f = open(snpsTxtFile)
	lines = f.readlines()
	f.close()
	
	currChr = '-1'
	strRes = ''
	for l in lines :
		if l[0] != '#' : #ignore comments
			sl = l.replace('\t\t', '\t').split('\t')
			if sl[0] != currChr :
				if currChr != '-1' :
					f = open('%s/%s.casavasnps'%(path, currChr), 'w')
					f.write(strRes)
					f.close()
				print 'importing snp data for %s...' % sl[0]
				currChr = sl[0]
				strRes = ''
			del(sl[0]) #remove chr
			del(sl[1]) #remove first position of the range
			strRes += ';'.join(sl)

	if currChr != '-1' :
		f = open('%s/%s.casavasnps'%(path, currChr), 'w')
		f.write(strRes)
		f.close()
	
	#print '%s/sourceFile.txt'%(path)
	f = open('%s/sourceFile.txt'%(path), 'w')
	f.write(snpsTxtFile)
	f.close()
	
	print 'importation of genome %s/%s done.' %(specie, genomeName)


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

	legend = ['//pos', 'chromosome', 'rs', 'type', 'alleles', 'validated', 'assembly', 'original_orientation', 'maf_allele', 'maf_count', 'maf', 'het', 'se(het)', 'loc']
	#Fcts
	def fillCSV(res, csv) :
		print '\t\t sorting snps by position in chromosome...'
		l = res.keys()
		l.sort()
		print '\t\t filling CSV...'
		for pos in l :
			li = csv.addLine()
			for field in csv.legend :
				csv.set(li, field, res[pos][field.lower()])
		
	def parse(s, chroNumber, legend) :
		lines = s.split('\n')
		
		res = {}
		
		for field in legend :
			res[field] = None
			
		criticalFields = ['rs', 'chromosome', '//pos', 'alleles', 'assembly', 'validated']
		numericFields = ['maf_count', 'maf', 'het', 'se(het)']
		numericFieldsWithNonNumericValues = 0
		locs = {}
		for l in lines :
			sl = l.split('|')
			if sl[0][:2] == 'rs' :
				res['rs'] = sl[0][2:].strip()
				res['type'] = sl[3].strip()
			
			elif sl[0][:3] == 'SNP' and res['rs'] != None :
				res['alleles'] = sl[1].strip().replace('alleles=', '').replace("'", "")
				res['het'] = sl[2].strip().replace('het=', '')				
				res['se(het)'] = sl[3].strip().replace('se(het)=', '')
				
			elif sl[0][:3] == 'VAL' and res['rs'] != None :
				res['validated'] = sl[1].strip().replace("validated=", '')
				
			elif sl[0][:3] == 'CTG' and sl[1].find('GRCh') > -1 and res['rs'] != None and (res['chromosome'] == None or res['//pos'] == None):
				res['original_orientation'] = sl[-1].replace('orient=', '').strip()
				res['assembly'] = sl[1].replace('assembly=', '').strip()
				res['chromosome'] = sl[2].replace('chr=', '').strip()
				pos = sl[3].replace('chr-pos=', '').strip()
				
				try:
					pos = int(pos)
					posOk = True
				except :
					posOk = False
					
				if res['chromosome'] != chroNumber or not posOk :
					res['chromosome'] = None
					res['//pos'] = None
					return (None, numericFieldsWithNonNumericValues)
				
				res['//pos'] = int(pos) -1
				
			elif sl[0][:4] == 'GMAF' and res['rs'] != None :
				res['maf_allele'] = sl[1].strip().replace('allele=', '')
				res['maf_count'] = sl[2].strip().replace('count=', '')
				res['maf'] = sl[3].strip().replace('MAF=', '')
			
			elif sl[0][:3] == 'LOC' and res['rs'] != None :
				allele = sl[4].strip().replace('allele=', '')
				locs[allele] = {}
				locs[allele]['gene'] = sl[1].strip()
				locs[allele]['function'] = sl[3].strip().replace('fxn-class=', '')
				try :
					locs[allele]['residue'] = sl[6].strip().replace('residue=', '')
				except :
					locs[allele]['residue'] = None
					
		for field in criticalFields :
			if res[field] == None :
				return (None, numericFieldsWithNonNumericValues)
		
		for field in numericFields :
			try :
				res[field] = float(res[field])
			except :
				res[field] = 0.0
				numericFieldsWithNonNumericValues += 1
				
		if res['original_orientation'] == '-' :
			res['alleles'] = uf.complement(res['alleles'])
		
		#res['loc'] = zlib.compress(pickle.dumps(locs))
		res['loc'] = pickle.dumps(locs).replace('\n', '/rje3/').replace(';', '/qte3/')
		return (res, numericFieldsWithNonNumericValues)
		
	#Fcts
	desc = 	"""pyGeno snp format is somewhat different from dbSNP's :
		- all SNPs that have a chromosome different from the one of the file, or chromosome number of '?' were discarded
		- the default value for numeric values (ex: het, se(het), 'MAF', 'GMAf count') is 0.0
		- no '?' allowed. If a numeric had a value of '?' (ex: het, se(het), 'MAF', 'GMAf count') the value was set to 0.0
		- all snps have an orientation of +. If a snp had an orientation of -, it's alleles were replaced by their complements
		- positions are 0 based
		- the only values extracted dbSNP files were from the files were : 'posistion', 'rs', 'type', 'assembly', 'chromosome', 'validated', 'alleles', 'original_orientation', 'maf_allele', 'maf_count', 'maf', 'het', 'se(het)'"""


	files = glob.glob(packageFolder+'/*.flat.gz')
	outPath = conf.pyGeno_SETTINGS['DATA_PATH']+'/%s/dbSNP/%s/' %(specie, versionName)
	if not os.path.exists(outPath):
		os.makedirs(outPath)

	for fil in files :
		chrStrStartPos = fil.find('ch')
		chrStrStopPos = fil.find('.flat')
		numericFieldsWithNonNumericValues = 0
		
		chroNumber = fil[chrStrStartPos+2: chrStrStopPos]
		outFile = fil.replace(packageFolder, outPath).replace('ds_flat_', '').replace('ch'+chroNumber, 'chr'+chroNumber).replace('.flat.gz', '.pygeno-dbSNP')
		
		print "extracting file :", fil, "..."
		f = gzip.open(fil)
		snps = f.read().split('\n\n')
		f.close()
		print snps[0]
		
		print "\tparsing..."
		
		resCSV = CSVFile(legend)
		res = {}
		for snp in snps[1:] :
			snpVals = parse(snp, chroNumber, legend)
			if snpVals[0] != None :
				res[snpVals[0]['//pos']] = snpVals[0]
			
			numericFieldsWithNonNumericValues += snpVals[1]
			
		print "\tformating data..."
		fillCSV(res, resCSV)

		print "\tsaving..."
		header = "source file : %s\n%s\n# of numeric fields with non mumeric values: %s (their values has been set to default 0.0)" % (fil, snps[0], numericFieldsWithNonNumericValues)
		header = header.split('\n')
		for i in range(len(header)) :
			header[i] = '//'+header[i]+'\n'
		header = ''.join(header)
		
		resCSV.setHeader(header)
		resCSV.save(outFile)

	readMeFile = outPath + 'README_format-description.txt'
	f = open(readMeFile, 'w')
	f.write(desc)
	f.close()

	
if __name__ == "__main__" :
	#importGenome('/u/daoudat/py/pyGeno/importationPackages/genomes/mouse', verbose = False)
	g = Genome(specie = 'Mus_musculus', name = 'GRCm38_test')
	print g.get(Exon, {'x1 >' : 797276})
	
	#print g.chromosomes
	#f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, Exon)
	#f.addFilter(**{'genome' : g, 'x1 <' : 797276})
	#for e in f.run() :
	#	print e.x1, e.x1 - 797276
