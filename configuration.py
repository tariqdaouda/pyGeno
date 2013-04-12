import sys, os, time

DATA_PATH = os.path.dirname(__file__)+ "/pyGenoData/"

pyGeno_BRANCH = "Medusa"
pyGeno_VERSION_NAME = 'Skiny Cobra'
pyGeno_VERSION_RELEASE_LEVEL = 'Beta'
pyGeno_VERSION_NUMBER = 13.04
pyGeno_VERSION_BUILD_TIME = time.ctime(os.path.getmtime(__file__))

pyGeno_VERSION_TUPLE = (pyGeno_BRANCH, pyGeno_VERSION_NAME, pyGeno_VERSION_RELEASE_LEVEL, pyGeno_VERSION_NUMBER, pyGeno_VERSION_BUILD_TIME )
pyGeno_VERSION_STR = "Branch: %s, Name: %s, Release Level: %s, Version: %s, Build time: %s" % pyGeno_VERSION_TUPLE

def configure():
	if sys.version_info[0] < 2 or (sys.version_info[0] > 2  and sys.version_info[1] < 7) :
		print "==> FATAL : pyGeno only works with python 2.7 and above, please upgrade your version"
		sys.exit(1)

configure()