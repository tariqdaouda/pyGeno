import os, glob, gzip, pickle, shutil, zlib
import configuration as conf

from Genome import Genome
from Chromosome import Chromosome
from Gene import Gene
from Transcript import Transcript
from Exon import Exon
from Protein import Protein

from tools import UsefulFunctions as uf
from tools.CSVTools import CSVFile

from expyutils.GTFTools import GTFFile

#from expyutils.GTFTools import GTFFile
	
"""===These are the only function that you will need to import new features  in pyGeno==="""
def importGenome(packageDir, specie, genomeName) :
	"""To import a new genome from scratch simply put all the compressed (in fa.gz) fasta files and the gtf (.gtf.gz) gene annotation file in the same folder. Then pass
	the path to this folder as the packageDir. Such files can be downloaded for example from ensembl: http://useast.ensembl.org/info/data/ftp/index.html
	
	In addition to thoses files the package should also contain a file named manifest.txt that contains the description of which file is associated to which chromosome,
	and the name of the gtf file.
	
	The format of the manifest is simply:
	>chormosomes
	1[\t]fileName
	2[\t]fileName
	3[\t]fileName
	....
	>gtf
	filename
	
	chromosomes must be in ordered in an inscreasing fashion
	You can add comments to the file with #. The manifest must end with an empty line
	
	If no manifest.txt is found in the package, pyGeno will try to create one. In order for this to work there should be only one .gtf.gz file in the package, and each chromosome
	fasta files (.fa.gz) must have the exact same nomenclature with only the chromosome name changing among the file names.
	"""
	
	try :
		f = open(packageDir+"/manifest.txt")
	except IOError:
		print "couldn't find the file manifest.txt in package"
		print "Here's the manifset:"
		print __createManifest(packageDir)
		
		while True :
			userInput = raw_input("\nIf the manifest is OK please enter 'yes'. If not please modify the file and then enter 'yes':\n")
			if userInput.upper() == 'YES' :
				break
	
	f = open(packageDir+"/manifest.txt")
	manifestLines = f.readlines()
	f.close()
	
	targetDir = conf.pyGeno_SETTINGS['DATA_PATH']+'/%s/genomes/%s'%(specie, genomeName)
	if not os.path.exists(targetDir):
		os.makedirs(targetDir)
	
	fastaHeaders = []
	chromosomePosIndex = []
	currGenomePos = 0
	
	parsing = None
	for line in manifestLines:
		if line == '\n' or line[0] == '#' :
			parsing = None
		elif line[:-1] == '>chromosomes' :
			parsing = 'chro'
		elif line[:-1] == '>gtf' :
			parsing = 'gtf'
			
		elif parsing == 'chro' :
			sl = line.split('\t')
			chroNumber = sl[0]
			res = __importSequence(chroNumber, "%s/%s" % (packageDir, sl[1][:-1]), targetDir)
			fastaHeaders.append(res[0])
			chromosomePosIndex.append("%s;%s;%s;%s" % (chroNumber, currGenomePos, currGenomePos + res[1], res[1]))
			currGenomePos = currGenomePos + res[1]
			
		elif parsing == 'gtf' :
			__importGeneSymbolIndex("%s/%s" % (packageDir, line[:-1]), targetDir)
	
	fastaHeadersFN = targetDir+'/fasta_headers.txt'
	print "saving fasta headers to %s..." % fastaHeadersFN
	f = open(fastaHeadersFN, 'w')
	f.write(''.join(fastaHeaders))
	f.close()
	
	genomeChroIndexFN = targetDir+'/genomeChrPos.index'
	print "saving saving chromosome position index to %s..." % genomeChroIndexFN
	f = open(genomeChroIndexFN, 'w')
	f.write('chromosome;chr_start;chr_end;length\n' + '\n'.join(chromosomePosIndex))
	f.close()	
	
	print "copying manifest to target directory..."
	shutil.copyfile(packageDir+"/manifest.txt", targetDir+"/manifest.txt")
	print 'done.'
>>>>>>> medusa:Importation.py


def __createManifest(packageDir) :
	
	fastas = glob.glob(packageDir+'/*.fa.gz')
	fastaStr = ">chromosomes\n"
	if len(fastas) == 1 :
		sf = fastas[0].replace(".fa.gz", '')
		number = sf[-1].upper()
		fastaStr += "%s\t%s" % (number, os.path.basename(fastas[0]))
	elif len(fastas) > 1 :
		commonSup = ""
		commonInf = ""
		i = 0
		add = True
		while add and i < len(fastas[0]) :
			for fasta in fastas[1:] :
				if i >= len(fasta) or fastas[0][i] != fasta[i] :
					add = False
					break
			if add :
				commonSup += fastas[0][i]
			i += 1
		
		i = 1
		add = True
		while add and i < len(fastas[0]) :
			for fasta in fastas[1:] :
				#print fastas[0][len(fastas[0])-i], fasta[len(fasta)-i], len(fastas[0])-i, len(fasta)-i
				if i < 0 or fastas[0][len(fastas[0])-i] != fasta[len(fasta)-i] :
					add = False
					break
			if add :
				commonInf += fastas[0][len(fastas[0])-i]
			i += 1
		commonInf = commonInf[::-1]
	
		chros = {}
		for fasta in fastas :
			number = fasta.replace(commonSup, '').replace(commonInf, '')
			try :
				chros[int(number)] = os.path.basename(fasta)
			except :
				chros[number] = os.path.basename(fasta)
			
		sortedNumbers = chros.keys()
		sortedNumbers.sort()
		fastaStr = ">chromosomes\n"
		for k in sortedNumbers :
			fastaStr += "%s\t%s\n" % (k, chros[k])

	gtfs = glob.glob(packageDir+'/*.gtf.gz')
	
	if len(gtfs) != 1 :
		raise Exception('There should be one and only one gtf index file in the package')
	
<<<<<<< HEAD:InstallationTools.py
	genome = Genome(name = genomeName, specie = specie, packageDir = packageDir)
	chroInfos = _importSequences(packageDir, genome)
	_importGenomeObjects(gtfs[0], chroInfos, genome)
	
	for chromosome in genomes.chromosomes :
		chromosome.header = chroInfos[chromosome.number][0]
		chromosome.x1 = chroInfos[chromosome.number][1]
		chromosome.x2 = chroInfos[chromosome.number][2]

	genome.save()
=======
	gtfStr = ">gtf\n%s" % os.path.basename(gtfs[0])
	
	comments = """#This file has been generated automaticly
#A manifest file must end with an empty line
#>[category] indicates the start of a new category of files to be parsed
#Empty lines are ok everywhere but within categories"""
	
	manifestStr = "%s\n%s\n\n%s\n" % (comments, fastaStr, gtfStr)
	f = open(packageDir+"/manifest.txt", "w")
	f.write(manifestStr)
	f.close()

	return manifestStr
>>>>>>> medusa:Importation.py
	
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
<<<<<<< HEAD:InstallationTools.py

"""===These are the private functions that you should not call, unless you really know what you're doing==="""
def _importSequences(fastaDir, genome) :
	print r"""Converting fastas from dir: %s into pyGeno's data format
	resulting files will be part of genome: %s/%s
	This may take some time, please wait...""" %(fastaDir, genome.specie, genome.name)
	

	path = conf.DATA_PATH+'/%s/genomes/%s/'%(genome.specie, genome.name)

	if not os.path.exists(path):
		os.makedirs(path)
		
	chrsFiles = glob.glob(fastaDir+'/chr*.fa')
	print chrsFiles, fastaDir+'/chr*.fa'
	if len(chrsFiles) < 1 :
		raise Exception('No Fasta found in directory, importation aborted')
		
	startPos = 0
	chroInfos = {} 
	headers = ''
	
	for chroFile in chrsFiles :
		print 'making sequence file for Chr', chroFile
		header = f.readline()
		
		f = open(chroFile)
		headers += f.readline()
		s = f.read().upper().replace('\n', '').replace('\r', '')
		f.close()
		
		fn = chroFile.replace('.fa', '.dat').replace(fastaDir, path)
		f = open(fn, 'w')
		f.write(s)
		f.close()
		
		startPos = startPos+len(s)
		chroInfos[chroNumber] = [header, str(startPos), str(startPos+len(s))]
	
	print 'writing build info file'
	f=open('%s/build_infos.txt' % (path), 'w')
	f.write('source package directory: %s' % fastaDir)
	f.write('\nheaders:\n------\n %s' % headers)
=======
	
def __importSequence(number, fastaFile, targetDir) :
	print 'making data for chromsome %s, source file: %s...' %(number, fastaFile)
	
	f = gzip.open(fastaFile)
	header = f.readline()
	strRes = f.read().upper().replace('\n', '').replace('\r', '')
>>>>>>> medusa:Importation.py
	f.close()
	
	fn = '%s/chr%s.dat' % (targetDir, number)#fastaFile.replace('.fa.gz', '.dat').replace(os.path.dirname(fastaFile), targetDir)
	f = open(fn, 'w')
	f.write(strRes)
	f.close()
	
	return (header, len(strRes))

<<<<<<< HEAD:InstallationTools.py
	return chroInfos
	
def _importGeneSymbolIndex_purgatory(gtfFile, specie) :

	path = conf.DATA_PATH+'/%s/gene_sets/'%(specie)
	if not os.path.exists(path):
		os.makedirs(path)
		
	f = open(gtfFile)
=======
def __importGeneSymbolIndex(gtfFile, targetDirectory) :
	
	f = gzip.open(gtfFile)
>>>>>>> medusa:Importation.py
	gtf = f.readlines()
	f.close()
	
	targetDir = targetDirectory+'/gene_sets/'
	
	if not os.path.exists(targetDir):
		os.makedirs(targetDir)
		
	chroLine = 0
	chro = -1
	currChro = -1
	gtfStr = ''
	pickIndex = {}
	symbol = -1
	currSymbol = -1
	for i in range(len(gtf)) :
		sl = gtf[i].split('\t')
		currChro = sl[0]
		if chro != currChro :
			if chro != -1 :
				print '\tsaving chromosme gtf file...'
				f = open('%s/chr%s.gtf' % (targetDir, chro), 'w')
				f.write(gtfStr)
				f.close()
				print '\tsaving chromosme pickled index...'
				f=open('%s/chr%s_gene_symbols.index.pickle' % (targetDir, chro), 'w')
				pickle.dump(pickIndex, f, 2)
				f.close()
				print '\tdone.'
			print 'making gene symbol index for chr %s...' %chro
			
			gtfStr = ''
			pickIndex = {}
			chro = currChro
			currSymbol = -1
			symbol = -1
			chroLine = 0
			
		currSymbol = gtf[i].split('\t')[8].split(';')[3]
		currSymbol = currSymbol[currSymbol.find('"') + 1:].replace('"', '').strip()
		
		if currSymbol != symbol :
			if symbol != -1 :
				pickIndex[symbol] = "%d;%d" %(startL, chroLine)
			
			symbol = currSymbol
			startL = chroLine
		
		gtfStr += gtf[i]
		chroLine += 1

	pickIndex[symbol] = "%d;%d" %(startL, chroLine)
	print '\tsaving chromosme gtf file...'
	f = open('%s/chr%s.gtf' % (targetDir, chro), 'w')
	f.write(gtfStr)
	f.close()
	print '\tsaving chromosme pickled index...'
	f=open('%s/chr%s_gene_symbols.index.pickle' % (targetDir, chro), 'w')
	pickle.dump(pickIndex, f, 2)
	f.close()
	print '\tdone.'
	

def _importGenomeObjects(gtfFilePath, chroInfos, genome) :

	gtf = GTFFile()
	gtf.parseFile(gtfFilePath)
	
	chromosomes = {}
	genes = {}
	transcripts = {}
	proteins = {}
	exons = {}
	for i in range(len(gtf)) :
		chroNumber = gtf.get(i, 'source')
		if chroNumber in chroInfos :
			if chroNumber not in chromosomes :
				print chroNumber, gtf.get(i, 'source')
				chromosomes[chroNumber] = Chromosome(genome = genome, number = chroNumber)
			
			geneId = gtf.get(i, 'gene_id')
			geneName = gtf.get(i, 'gene_name')
			strand = gtf.get(i, 'strand')
			if geneId not in genes :
				genes[geneId] = Gene(genome = genome, chromosome = chromosomes[chroNumber], id = geneId, name = geneName, strand = strand)
			
			transId = gtf.get(i, 'transcript_id')
			transName = gtf.get(i, 'transcript_name')
			protId = gtf.get(i, 'protein_id')
			protName = transName
			if transId not in transcripts :
				transcripts[transId] = Transcript(genome = genome, chromosome = chromosomes[chroNumber], gene = genes[geneId], id = transId, name = transName)
				protein[protId] = Protein(genome = genome, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], id = protId, name = protName)
				transcripts[transId].protein = protein[protId]
				
			regionType = gtf.get(i, 'feature')
			x1 = gtf.get(i, 'start')
			x2 = gtf.get(i, 'end')
			if regionType == 'exon' :
				exonNumber = gtf.get(i, 'exon_number')
				codingType = gtf.get(i, 'gene_biotype')
				exonKey = '%s%s%s' % (chroNumber, x1, x2)
				exonId = gtf.get(i, 'exon_id')
				if exonId not in exons :
					exons[exonKey] = Exon(genome = genome, chromosome = chromosomes[chroNumber], gene = genes[geneId], transcript = transcripts[transId], strand = strand, number = exonNumber, gene_biotype = gene_biotype, id = exonId)
			
			elif regionType == 'CDS' :
				exons[exonKey].setCDS(x1, x2)
			elif regionType == 'start_codon' :
				exons[exonKey].startCodon = x1
			elif regionType == 'stop_codon' :
				exons[exonKey].stopCodon = x2
		
		for transcript in transcripts :		
			f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, Exon)
			f.addFilter(**{'transcript' : transcript})
			transcript.exons = f.run()
		
		for gene in genes :
			f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, Transcript)
			f.addFilter(**{'gene' : gene})
			gene.transcripts = f.run()
			
			f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, Exon)
			f.addFilter(**{'gene' : gene})
			gene.exons = f.run()
			
		for chro in chromosomes :
			f = RabaQuery(conf.pyGeno_RABA_NAMESPACE, Gene)
			f.addFilter(**{'chromosome' : chro})
			chro.genes = f.run()
		
	genome.chromosomes = RabaList(chromosomes.values())
	genome.sourceGTF = gtfFilePath
	
if __name__ == "__main__" :
<<<<<<< HEAD:InstallationTools.py
	#import_dbSNP('/u/daoudat/py/pyGeno/importationPackages/dbSNP/human/dbSNP137', 'human', 'dbSNP137')
=======
	#importGenome("/u/daoudat/py/pyGeno/installationPackages/genomes/mouse", 'mouse', 'GRCm38')
	import_dbSNP('/u/daoudat/py/pyGeno/importationPackages/dbSNP/human/dbSNP137', 'human', 'dbSNP137')
>>>>>>> medusa:Importation.py
	#importGenome_casava('human', 'lightR_Transcriptome', '/u/corona/Project_DSP008a/Build_Diana_ARN_R/snps.txt')
	print 'import mouse'
	importGenome('/u/daoudat/py/pyGeno/importationPackages/genomes/mouse', 'mouse', 'mouseRaba')
	#print 'import b6'
	#importGenome_casava('mouse', 'B6', '/u/corona/Project_DSP014/120313_SN942_0105_AD093KACXX/Build_B6/snps.txt')
	
