import os, glob, gzip, pickle
import configuration as conf

#import cPickle
from Genome import Genome
from tools import UsefulFunctions as uf
from expyutils.CSVTools import CSVFile
#from expyutils.GTFTools import GTFFile

def installSequences(fastaDir, specie, genomeName) :
	print r"""Converting fastas from dir: %s into pyGeno's data format
	resulting files will be part of genome: %s/%s
	This may take some time, please wait...""" %(fastaDir, specie, genomeName)
	
	
	#path = conf.DATA_PATH+'/ncbi/%s/sequences/%s/'%(specie, genomeName)
	path = conf.DATA_PATH+'/%s/genomes/%s/'%(specie, genomeName)

	if not os.path.exists(path):
		os.makedirs(path)
		
	chrs = glob.glob(fastaDir+'/chr*.fa')
	if len(chrs) < 1 :
		raise Exception('No Fastas found in directort, installation aborted')
		
	#chrs.extend(glob.glob(DATA_PATH/+'*.fasta'))
	startPos = 0
	genomeChrPos = {} 
	headers = ''
	
	for chro in chrs :
		print 'making data for Chr', chro
		
		f = open(chro)
		headers += f.readline()
		s = f.read().upper().replace('\n', '').replace('\r', '')
		f.close()
		
		fn = chro.replace('.fa', '.dat').replace(fastaDir, path)
		f = open(fn, 'w')
		f.write(s)
		f.close()
		
		ch = chro.replace(fastaDir, '').replace('/chr', '').replace('.fa', '')
		genomeChrPos[ch] = [str(ch), str(startPos), str(startPos+len(s)), str(len(s))]
		startPos = startPos+len(s)
		
	#fh = open(path+'/fasta_headers.txt', 'w')
	#fh.write(headers)
	#fh.close()

	print 'making index file genomeChrPos'
	fh = open(path+'/genomeChrPos.index', 'w')
	fh.write('chromosome;chr_start;chr_end;length\n')
	for k in genomeChrPos.keys() :
		fh.write(';'.join(genomeChrPos[k])+'\n')
	fh.close()

	print 'writing build info file'
	f=open('%s/build_infos.txt' % (path), 'w')
	f.write('source package directory: %s' % fastaDir)
	f.write('\nheaders:\n------\n %s' % headers)
	f.close()

def installGeneSymbolIndex(gtfFile, specie) :
	#gtf = GTFFile()
	#gtf.parseFile(gtfFile)
	
	path = conf.DATA_PATH+'/%s/gene_sets/'%(specie)
	if not os.path.exists(path):
		os.makedirs(path)
		
	f = open(gtfFile)
	gtf = f.readlines()
	f.close()
	
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
				f = open('%s/chr%s.gtf' % (path, chro), 'w')
				f.write(gtfStr)
				f.close()
				print '\tsaving chromosme pickled index...'
				f=open('%s/chr%s_gene_symbols.index.pickle' % (path, chro), 'w')
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
		#print currSymbol
		
		if currSymbol != symbol :
			if symbol != -1 :
				pickIndex[symbol] = "%d;%d" %(startL, chroLine)
			
			symbol = currSymbol
			startL = chroLine
		
		gtfStr += gtf[i]
		chroLine += 1

	pickIndex[symbol] = "%d;%d" %(startL, chroLine)
	print '\tsaving chromosme gtf file...'
	f = open('%s/chr%s.gtf' % (path, chro), 'w')
	f.write(gtfStr)
	f.close()
	print '\tsaving chromosme pickled index...'
	f=open('%s/chr%s_gene_symbols.index.pickle' % (path, chro), 'w')
	pickle.dump(pickIndex, f, 2)
	f.close()
	print '\tdone.'
	print 'writing build info file'
	f=open('%s/build_infos.txt' % (path), 'w')
	f.write('source file: %s' % gtfFile)
	f.close()

def installGenome(packageDir, specie, genomeName) :
	gtfs = glob.glob(packageDir+'/*.gtf')
	print gtfs
	if len(gtfs) != 1 :
		raise Exception('There should be one and only one gtf index file in the package')
	
	installSequences(packageDir, specie, genomeName)
	installGeneSymbolIndex(gtfs[0], specie)
	
def installGenome_casava(specie, genomeName, snpsTxtFile) :
	"""Creates a light genome (contains only snps infos and no sequence from the reference genome)
	The .casavasnps files generated are identical to the casava snps but with ';' instead of tabs and 
	a single position instead of a range"""

	#path = conf.DATA_PATH+'/ncbi/%s/sequences/%s/'%(specie, genomeName)
	path = conf.DATA_PATH+'/%s/genomes/%s/'%(specie, genomeName)
	print 'Installing genome %s...' %path
	
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
				print 'Installing snp data for %s...' % sl[0]
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
	
	print 'Installation of genome %s/%s done.' %(specie, genomeName)

def makeDbSNPFile_chrrpts_bck(filePath, outputFolder, compressOutput = False, verbose = False) :
	"""
	filePath should be a gziped file.
	
	The resulting file is almost identical to dbSNPs except that:
	-they are ';' separated 
	-dbSNP Files are sorted by rsId, in pyGeno format they are sorted by position in the chromosome. This function does the sorting.
	-In pygeno files the position in the chormosome and rsId column are inverted => chr position is now at column 0 and the rsid at column 11
	-The chr prostion begins at 0 and not a 1 as in dbSNP
	
	If compressOutput is set to true, the sorted output for each chromosome will be a gz file.
	This creates significantly smaller files on disk but will result in an overhead while loading the data"""
	
	print "Extracting data from %s..." %filePath
	f = gzip.open(filePath, 'rb')#open('chr_Y.txt')
	lines = f.readlines()
	f.close()
	
	d = {}
	print 'populating dict...'
	for l in lines :
		sl = l.split('\t');
		try :
			chroPos = str(int(sl[11])-1) #pyGeno starts counting at 0
			sl[11] = sl[0]
			sl[0] = chroPos
			d[sl[11]] = ';'.join(sl)
		except IndexError:
			if verbose:
				print "ignoring invalid line %s" %l
		except ValueError:
			if verbose:
				print "ignoring line for invalid chromosome position : ", sl[11]
			
	print 'sorting...'
	keys = d.keys()
	keys.sort()
	
	filename = filePath.split('/')[-1]
	outFile = filename.replace('chr_', 'chr')
	outFile = outFile.replace('.txt.gz', '.pygeno-dbsnp')
	outFile = outputFolder + '/' + outFile
	if compressOutput :
		if verbose :
			print "\toutput will is compressed"
		print "writing output to %s..." % (outFile+'.gz')
		f = gzip.open(outFile+'.gz', 'w')
	else :
		print "writing output to %s..." % (outFile)
		f = open(outFile, 'w')

	for k in keys :
		f.write(d[k])
	f.close()

def install_dbSNP(packageFolder, specie, versionName) :
	"""To install dbSNP informations, download ASN1_flat files from the 
	dbSNP ftp : ftp://ftp.ncbi.nih.gov/snp/organisms/ and place them all in one single folder. This folder
	will be considered as a package.
	Launch this function and go make yourself a cup of coffee, this function has absolutly not been written to be fast
	
	versionName is name with wich you want to call this specific version of dbSNP
	"""
	#Fcts
	def fillCSV(res, csv) :
		l = res.keys()
		l.sort()
		for k in l :
			li = csv.addLine()
			csv.setElement(li, '//pos', res[k]['pos'])
			csv.setElement(li, 'chro', res[k]['chro'])
			csv.setElement(li, 'rs', res[k]['rs'])
			csv.setElement(li, 'type', res[k]['type'])
			csv.setElement(li, 'alleles', res[k]['alleles'])
			csv.setElement(li, 'validated', res[k]['validated'])
			csv.setElement(li, 'assembly', res[k]['assembly'])
		
	def parse(s, chroNumber, res) :
		lines = s.split('\n')
		
		rs = None
		typ = None
		pos = None
		assembly = None
		chro = None
		validated = None
		alleles = None
		
		for l in lines :
			sl = l.split('|')
			if sl[0][:2] == 'rs' :
				rs = sl[0][2:].strip()#csv.setElement(li, 'rs', )
				typ = sl[3].strip()#csv.setElement(li, 'rs', )
			
			elif sl[0][:3] == 'SNP' and rs != None :
				alleles = sl[1].strip().replace('alleles=', '').replace("'", "")
				#csv.setElement(li, 'alleles', )
			
			elif sl[0][:3] == 'VAL' and rs != None :
				validated = sl[1].strip().replace("validated=", '')
				#csv.setElement(li, 'alleles', )
				
			elif sl[0][:3] == 'CTG' and sl[1].find('GRCh') > -1 and rs != None :
				assembly = sl[1].replace('assembly=', '').strip()
				chro = sl[2].replace('chr=', '').strip()
				pos = sl[3].replace('chr-pos=', '').strip()
				if chro != chroNumber or pos == '?' :
					chro = None
					pos = None
				else :
					pos = int(pos)-1
				
			if rs != None and chro != None and pos != None and alleles != None and assembly != None and validated != None:
				res[pos] = {'pos' : pos, 'chro' : chro, 'rs' : rs, 'type' : typ, 'alleles' : alleles, 'validated' : validated, 'assembly' : assembly}
				break
	#Fcts
	files = glob.glob(packageFolder+'/*.flat.gz')
	outPath = conf.DATA_PATH+'/%s/dbSNP/%s/' %(specie, versionName)
	if not os.path.exists(outPath):
		os.makedirs(outPath)

	for fil in files :
		chrStrStartPos = fil.find('ch')
		chrStrStopPos = fil.find('.flat')
		#print chrStrStartPos, chrStrStopPos
		chroNumber = fil[chrStrStartPos+2: chrStrStopPos]
		outFile = fil.replace(packageFolder, outPath).replace('ds_flat_', '').replace('ch'+chroNumber, 'chr'+chroNumber).replace('.flat.gz', '.pygeno-dbSNP')
		headerFile = outFile.replace('.pygeno-dbSNP', '-header.txt')
		
		print "extracting file :", fil, "..."
		f = gzip.open(fil)
		snps = f.read().split('\n\n')
		f.close()
		print snps[0]
		print "\tparsing..."
		res = {}
		resCSV = CSVFile(['//pos', 'chro', 'rs', 'type', 'alleles', 'validated', 'assembly'])
		for snp in snps[1:] :
			parse(snp, chroNumber, res)
			
		print "\tformating data..."
		fillCSV(res, resCSV)

		print "\tsaving..."
		header = "source file : %s\n%s" % (fil, snps[0])
		header = header.split('\n')
		for i in range(len(header)) :
			header[i] = '//'+header[i]+'\n'
		header = ''.join(header)
		
		resCSV.setHeader(header)
		resCSV.save(outFile)
		
		
		f = open(headerFile, 'w')
		f.write(header)
		f.close()
	
if __name__ == "__main__" :
	#installGenome("~/py/mous", 'antoine', 'tariq')
	#install_dbSNP('/u/daoudat/py/pyGeno/pyGenoData/installationPackages/dbSNP/human/dbSNP137', 'human', 'dbSNP137')
	#installGenome_casava('human', 'lightR_Transcriptome', '/u/corona/Project_DSP008a/Build_Diana_ARN_R/snps.txt')
	#makeGenome_casava('human', 'lightM_Transcriptome', '/u/corona/Project_DSP008a/Build_Diana_ARN_M/snps.txt')
	#makeGenome_casava('human', 'lightR_Exome', '/u/corona/Project_DSP008a/Build_Diana_ADN_R/snps.with_removed.txt')
	#makeGenome_casava('human', 'lightM_Exome', '/u/corona/Project_DSP008a/Build_Diana_ADN_M/snps.with_removed.txt')

	#print 'install mouse'
	#installGenome('/u/daoudat/py/pyGeno/mouse', 'mouse', 'reference2')
	#print 'install b6'
	#makeGenome_casava('mouse', 'B6', '/u/corona/Project_DSP014/120313_SN942_0105_AD093KACXX/Build_B6/snps.txt')
	#makeGeneSymbolIndex('pyGenoData/installs/mouse/Mus_musculus.NCBIM37.64.gtf', 'mouse')
	#installGenome('/u/daoudat/py/pyGeno/pyGenoData/installationPackages/mouse', 'mouse', 'reference')

