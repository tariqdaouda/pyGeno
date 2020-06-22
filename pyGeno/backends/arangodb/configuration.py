from ..backend_abs import DatabaseConfiguration_ABS
import pyArango.connection as CON
import pyArango.theExceptions as PEXP 

from . import savers
from . import query_handler

# _DB_CONNECTION = None
# _DB = None

class DatabaseConf(DatabaseConfiguration_ABS):
    """
    Configuration for using ArangoDB
    """
    _DB_NAME = "pyGeno"

    def reset(self, pyarango_conn_args) :
        """This function is automatically called at launch"""
        self.conn_args = pyarango_conn_args
        self.connection = CON.Connection(**self.conn_args)
        try :
            self.database = self.connection[self._DB_NAME]
        except :
            self.set_database()
        self.database = self.connection[self._DB_NAME]

    def list_genomes(self):
        """"""
        aql = """
            FOR g IN Genome
                RETURN DISTINCT g.name
        """
        ret = self.database.AQLQuery(aql, batchSize=1000)
        return ret.result

    def get_query_handler(self):
        self.query_handler = query_handler.QueryHandler(self)
        return self.query_handler
            
    def load_saver(self):
        self.saver = savers.GenomeSaver(self)
        return self.saver
    
    def set_database(self):
        try :
            self.connection.createDatabase(self._DB_NAME)
        except Exception as e :
            print("Unable to create database: %s" % e)
    
    def create_collection(self, colname):
        try :
            self.db.createCollection(colname)
        except PEXP.CreationError as e :
            print("Unable to create collection '%s' : %s" % (colname, e))

    def get_link_collection_name(self, col1, col2):
        if col1[0] > col2[0] :
            return "%s_X_%s" % (col1, col2)
        return "%s_X_%s" % (col2, col1)

    def get_configuration(self):
        """return configuration of the backend in a dictionary"""
        return {
            "module": "pyGeno.backends.arangodb",
            "arguments": self.conn_args
        }

    def prompt_setup(self):
        args = {
            "username": "root",
            "password": "root",
            "arangoURL": "http://localhost:8529/"
        }
        print("## ArangoDB backend setup")
        for key, value in args.items():
            val = input("~:> %s? (Press enter to use default is '%s'): " % (key, value))
            if val != "":
                args[key] = val

        self.reset(args)
