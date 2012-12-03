import os, glob, pickle
import configuration as conf
#import cPickle

from Genome import Genome
from tools import UsefulFunctions as uf
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
	
def makeGenome_casava(specie, genomeName, snpsTxtFile) :
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

def makeDbSNPFile(filePath, outputFolder, compressOutput = False, verbose = False) :
	"""
	filePath should be a gziped file.
	dbSNP Files are sorted by rsId, in pyGeno format they are sorted by position in the chromosome.
	This function does the sorting.
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
			sl[11] = str(int(sl[11])-1) #pyGeno starts counting at 0
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

def install_dbSNP(packageFolder, specie, versionName, compressOutput = False, verbose = False) :
	"""To install dbSNP informations, download chromosome reports (chr_rpts )files from the 
	dbSNP ftp : ftp://ftp.ncbi.nih.gov/snp/organisms/ and place them all in one single folder. This folder
	will be considered as a package. Launch this function and drink a cup of coffee while waiting
	
	versionName is name with wich you want to call this specific version of dbSNP
	If compressOutput is set to true, the sorted output for each chromosome will be a gz file.
	This creates significantly smaller files on disk but will result in an overhead while loading the data
	"""
	
	files = glob.glob(packageFolder+'/*.txt.gz')
	outputPath = conf.DATA_PATH+'/%s/dbSNP/%s/' %(specie, versionName)
	
	if not os.path.exists(outputPath):
		os.makedirs(outputPath)
		
	for fil in files :
		makeDbSNPFile(fil, outputPath)
	
if __name__ = "__main__" :
	#makeGenome_casava('human', 'lightR_Transcriptome', '/u/corona/Project_DSP008a/Build_Diana_ARN_R/snps.txt')
	#makeGenome_casava('human', 'lightM_Transcriptome', '/u/corona/Project_DSP008a/Build_Diana_ARN_M/snps.txt')
	#makeGenome_casava('human', 'lightR_Exome', '/u/corona/Project_DSP008a/Build_Diana_ADN_R/snps.with_removed.txt')
	#makeGenome_casava('human', 'lightM_Exome', '/u/corona/Project_DSP008a/Build_Diana_ADN_M/snps.with_removed.txt')

	#print 'install mouse'
	#installGenome('/u/daoudat/py/pyGeno/pyGenoData/fasta-installs/mouse', 'mouse', 'reference')
	#print 'install b6'
	#makeGenome_casava('mouse', 'B6', '/u/corona/Project_DSP014/120313_SN942_0105_AD093KACXX/Build_B6/snps.txt')
	#makeGeneSymbolIndex('pyGenoData/installs/mouse/Mus_musculus.NCBIM37.64.gtf', 'mouse')
	#installGenome('/u/daoudat/py/pyGeno/pyGenoData/installationPackages/mouse', 'mouse', 'reference')

