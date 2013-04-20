import sys, os, time, cPickle
from tools import io

pyGeno_SETTINGS_DIR = '~/.pyGeno/'
pyGeno_SETTINGS_DIR = os.path.expanduser(pyGeno_SETTINGS_DIR)
pyGeno_SETTINGS_FILENAME = '%s/pyGeno_settings.pkl' % pyGeno_SETTINGS_DIR

pyGeno_BRANCH = "Medusa"
pyGeno_VERSION_NAME = 'Skiny Cobra'
pyGeno_VERSION_RELEASE_LEVEL = 'Beta'
pyGeno_VERSION_NUMBER = 13.04
pyGeno_VERSION_BUILD_TIME = time.ctime(os.path.getmtime(__file__))

def version(self) :
	"""returns a tuple describing pyGeno's current version"""
	return (pyGeno_BRANCH, pyGeno_VERSION_NAME, pyGeno_VERSION_RELEASE_LEVEL, pyGeno_VERSION_NUMBER, pyGeno_VERSION_BUILD_TIME )

def prettyVersion(self) :
	"""returns pyGeno's current version in a pretty human redable way"""
	return "---\npyGeno ~-~-~-:> Branch: %s, Name: %s, Release Level: %s, Version: %s, Build time: %s" % pyGeno_VERSION_TUPLE, '\n---'

def checkPythonVersion(self) :
	if sys.version_info[0] < 2 or (sys.version_info[0] > 2  and sys.version_info[1] < 7) :
		return False
	return True

def checkDataPath(self) :
	if pyGeno_SETTINGS['DATA_PATH'] == None or not os.path.exists(pyGeno_SETTINGS['DATA_PATH']) :
		return False
	return True

def checkReferenceGenome(self, specie) :
	if specie not in pyGeno_SETTINGS['REFERENCE_GENOMES'] or not os.path.exists(pyGeno_SETTINGS['DATA_PATH']+'/genomes/'+pyGeno_SETTINGS['REFERENCE_GENOMES'][specie]) :
		return False
	return True

def change_REFERENCE_GENOME(self, specie) :
	strGenomes = ''
	genomeNames = set()
	try :
		for d in os.walk(pyGeno_SETTINGS['DATA_PATH']+'/%s/genomes/' % specie).next()[1] :
			strGenomes += '\t' + d + '\n'
			genomeNames.add(d)
	except StopIteration :
		print "There's no available genomes for specie: %s. Please import one using ImportTools.importGenome.\n" % specie
		return False
	
	stop = False
	while not stop :
		enterMsg = "Please enter the name of the genome you want to use as reference for %s, the possible choices are:\n%s\nYour choice: " % (specie, strGenomes)
		ret = io.enterConfirm(enterMsg)			
		if ret == None :
			return False
		
		if ret in genomeNames :
			stop = True
		else :
			print "%s is not available, please enter a valid name." % ret
	
	pyGeno_SETTINGS['REFERENCE_GENOMES'][specie] = ret
	cPickle.dump(pyGeno_SETTINGS, open(pyGeno_SETTINGS_FILENAME, 'w'))
	
	return True
	
def change_DATA_PATH(self) :
	enterMsg = "Please enter the absolute path you want to store the databases into (if it doesn't exist i'll try to create it for you).\nYou can always change it using change_DATA_PATH() :\n"
	ret = io.enterConfirm(enterMsg)
	if ret == None :
		sys.exit(1)
	pyGeno_SETTINGS['DATA_PATH'] = ret

	if not os.path.exists(pyGeno_SETTINGS['DATA_PATH'] ) :
		os.makedirs(pyGeno_SETTINGS['DATA_PATH'])
	pyGeno_SETTINGS['DATA_PATH'] = os.path.expanduser(pyGeno_SETTINGS['DATA_PATH'])
	cPickle.dump(pyGeno_SETTINGS, open(pyGeno_SETTINGS_FILENAME, 'w'))
	
def configure(self) :
	"This function is automaticly called at launch"
	
	if not checkPythonVersion() :
		print "==> FATAL : pyGeno only works with python 2.7 and above, please upgrade your version"
		sys.exit(1)

	if not os.path.exists(pyGeno_SETTINGS_DIR) :
		os.makedirs(os.path.expanduser(pyGeno_SETTINGS_DIR))
	
	try :
		pyGeno_SETTINGS = cPickle.load(open(pyGeno_SETTINGS_FILENAME))
	except :
		pyGeno_SETTINGS = {'DATA_PATH' : '', 'REFERENCE_GENOMES' : {}}
		cPickle.dump(pyGeno_SETTINGS, open(pyGeno_SETTINGS_FILENAME, 'w'))

	if not checkDataPath() :
		print "==> The pyGeno data folder is not set, or doesn't exist. The current value is: %s" % pyGeno_SETTINGS['DATA_PATH']
		change_DATA_PATH()

if __name__== '__main__' :
	change_REFERENCE_GENOME('human')#configure()

