import sys, os, time

# from configparser import SafeConfigParser

# class PythonVersionError(Exception) :
#     pass

_FACE = "~-~-:>"
_CHAPTER = "3"

_VERSION_NAME = 'Mighty Cobra'
_VERSION_RELEASE_LEVEL = 'Alpha'
_VERSION_NUMBER = 20.03
_VERSION_BUILD_TIME = time.ctime(os.path.getmtime(__file__))

_SETTINGS_DIR = os.path.normpath(os.path.expanduser('~/.pyGeno/'))
_DB_CONF_FILE = os.path.join(_SETTINGS_DIR, "db_conf.json")
_REMOTE_LOCATION = 'http://bioinfo.iric.ca/~feghalya/pyGeno_datawraps'

_BACKEND = None

def system_message(msg):
    sys.stderr.write(_FACE + " " + msg + "\n")

def version() :
    """returns a tuple describing pyGeno's current version"""
    return (_CHAPTER, _VERSION_NAME, _VERSION_RELEASE_LEVEL, _VERSION_NUMBER, _VERSION_BUILD_TIME )

def pretty_version() :
    """returns pyGeno's current version in a pretty human readable way"""
    return "pyGeno Chapter: %s, Name: %s, Release Level: %s, Version: %s, Build time: %s" % version()

def set_backend(backend=None, make_default=False):
    """
    Set the backend for pyGeno. if backend=None, pyGeno will revert to the
    default engine: RabaDB. If make_default, the backend will be stored
    as the new default.
    """
    from pyGeno.backends.rabadb.configuration import DatabaseConf
    import json

    global _BACKEND
    global _SETTINGS_DIR
    global _DB_CONF_FILE

    if backend is not None:
        sys.stderr.write("Switching to custom backend.")
        _BACKEND = backend
    else:
        if backend is None and _BACKEND is None:
            if os.path.exists(_DB_CONF_FILE):
                system_message("Using user defined default backend.")
                with open(_DB_CONF_FILE, 'r') as file:
                    json_conf = json.load(file)
                backend_module = importlib.import_module(json_conf["python_module"]+".configuration")
                _BACKEND = backend_module(**json_conf["arguments"])
            else:
                system_message("Using default RabaDB backend.")
                _BACKEND = DatabaseConf(_SETTINGS_DIR) 

    if make_default:
        with open(_DB_CONF_FILE, 'w') as file:
            json_conf = json.dump(_BACKEND.get_configuration(), file)
        system_message("Saved new default backend configuration.")

def get_backend():
	return _BACKEND

def pyGeno_init():
    global _SETTINGS_DIR
    if not os.path.isdir(_SETTINGS_DIR):
        os.mkdir(_SETTINGS_DIR)
    system_message(prettyVersion())
    set_backend()

if __name__ == '__main__':
	pyGeno_init