from ..backend_abs import DatabaseConfiguration_ABS
import pyArango.connection as CON

_DB_CONNECTION = None
_DB = None

class DatabaseConf(DatabaseConfiguration_ABS):
    """
    Configuration for using ArangoDB
    """
    _DB_NAME = "pyGeno"

    def reset(self, **pyarango_conn_args) :
        """This function is automatically called at launch"""
        self.conn_args = pyarango_conn_args
        self.connection = CON.Connection(**self.conn_args)
        self.database = self.connection[self._DB_NAME]

    def set_database(self):
        try :
          self.db = self.connection.createDatabase(self._DB_NAME)
        except Exception as e :
          print("Unable to create database: %s" % e)
          self.db = self.connection[self._DB_NAME]

        for colname in ("Peptides", "VirusSequences"):
          try :
            self.db.createCollection(colname)
          except PEXP.CreationError as e :
            print("Unable to create collection '%s' : %s" % (colname, e))

    def get_configuration(self):
        """return configuration of the backend in a dictionary"""
        return {
            "module": "pyGeno.backends.arangodb",
            "arguments": self.conn_args
        }