import sys, os, time
import rabaDB.setup

import user_config as uc

class PythonVersionError(Exception) :
	pass

pyGeno_FACE = "~-~-:>"
pyGeno_BRANCH = "V2"

pyGeno_VERSION_NAME = 'Lean Viper!'
pyGeno_VERSION_RELEASE_LEVEL = 'Beta'
pyGeno_VERSION_NUMBER = 14.02
pyGeno_VERSION_BUILD_TIME = time.ctime(os.path.getmtime(__file__))

pyGeno_RABA_NAMESPACE = 'pyGenoRaba'
pyGeno_RABA_DBFILE = os.path.normpath('%s/pyGenoRaba.db' % uc.pyGeno_SETTINGS_PATH)
pyGeno_DATA_PATH = os.path.normpath('%s/data' % uc.pyGeno_SETTINGS_PATH)

db = None

def version() :
	"""returns a tuple describing pyGeno's current version"""
	return (pyGeno_FACE, pyGeno_BRANCH, pyGeno_VERSION_NAME, pyGeno_VERSION_RELEASE_LEVEL, pyGeno_VERSION_NUMBER, pyGeno_VERSION_BUILD_TIME )

def prettyVersion() :
	"""returns pyGeno's current version in a pretty human redable way"""
	return "pyGeno %s Branch: %s, Name: %s, Release Level: %s, Version: %s, Build time: %s" % version()

def checkPythonVersion() :
	#pyGeno needs python 2.7+
	if sys.version_info[0] < 2 or (sys.version_info[0] > 2  and sys.version_info[1] < 7) :
		return False
	return True

def getGenomeSequencePath(specie, name) :
	return os.path.normpath(pyGeno_DATA_PATH+'/%s/%s' % (specie.lower(), name))

def pyGeno_init() :
	"This function is automaticly called at launch"
	
	global db
	
	if not checkPythonVersion() :
		raise PythonVersionError("==> FATAL: pyGeno only works with python 2.7 and above, please upgrade your python version")

	if not os.path.exists(uc.pyGeno_SETTINGS_PATH) :
		os.makedirs(uc.pyGeno_SETTINGS_PATH)

	if not os.path.exists(pyGeno_DATA_PATH) :
		os.makedirs(pyGeno_DATA_PATH)

	#launching the db
	rabaDB.setup.RabaConfiguration(pyGeno_RABA_NAMESPACE, pyGeno_RABA_DBFILE)
	db = rabaDB.setup.RabaConnection(pyGeno_RABA_NAMESPACE)
