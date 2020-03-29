import os
import rabaDB.rabaSetup
import rabaDB.Raba
from ..backend_abs import DatabaseConfiguration_ABS

class DatabaseConf(DatabaseConfiguration_ABS):
    """
    Configuration for using legacy RabaDB backend.
    Annotations are saved in a sqlite database and sequences
    are saved separately on disk and loaded using a memory map.
    This backend is slower, and handles less data but it doesn't
    require the installation of a third party database.
    """
    def reset(self, settings_path) :
        """This function is automatically called at launch"""
        
        self.NAMESPACE = 'pyGenoRaba'
        
        self.SETTINGS_PATH = settings_path
        self.DBFILE = os.path.normpath( os.path.join(self.SETTINGS_PATH, "pyGenoRaba.db") )
        self.DATA_PATH = os.path.normpath( os.path.join(self.SETTINGS_PATH, "data") )
        
        if not os.path.exists(self.DATA_PATH) :
            os.makedirs(self.DATA_PATH)

        #launching the db
        rabaDB.rabaSetup.RabaConfiguration(self.NAMESPACE, self.DBFILE)
        self.db = rabaDB.rabaSetup.RabaConnection(self.NAMESPACE)
        self.dbConf = rabaDB.rabaSetup.RabaConfiguration(self.NAMESPACE)

    def getGenomeSequencePath(self, specie, name) :
        return os.path.normpath(self.DATA_PATH+'/%s/%s' % (specie.lower(), name))

    def removeFromDBRegistery(self, obj) :
        """rabaDB keeps a record of loaded objects to ensure consistency between different queries.
        This function removes an object from the registery"""
        rabaDB.Raba.removeFromRegistery(obj)

    def freeDBRegistery(self, ) :
        """rabaDB keeps a record of loaded objects to ensure consistency between different queries. This function empties the registery"""
        rabaDB.Raba.freeRegistery()

    def get_configuration(self):
        """return configuration of the backend in a dictionary"""
        return {
            "module": "pyGeno.backends.rabadb",
            "arguments": {
                "settings_path": self.SETTINGS_PATH
            }
        }