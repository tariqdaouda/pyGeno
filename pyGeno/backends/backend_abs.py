class DatabaseConfiguration_ABS:
    """The general interface for objects storing database configurations"""
    def __init__(self, conf=None):
        if conf:
            self.reset(conf)
        self.saver = None

    def load_saver():
        raise NotImplemented("This is an abstract class")

    def reset(self, conf):
        raise NotImplemented("This is an abstract class")

    def get_configuration(self):
        """
        return configuration of the backend in a dictionary
        in the following format:

        {
            "module": "pyGeno.backends.rabadb",
            "arguments": {
                "settings_path": self.SETTINGS_PATH
            }
        }
        
        The return value can be stored by pyGeno in a json format
        to be used as default.
        """
        raise NotImplemented("This is an abstract class")
    
    def prompt_setup(self):
        raise NotImplemented("This is an abstract class")

class GenomeSaver_ABS(object):
    """
    Saves genome into database
    """
    def __init__(self, database_configuration):
        super(GenomeSaver_ABS, self).__init__()
        self.database_configuration = database_configuration
        self.store = {}
    
    def add(self, obj_type, unique_id, dct_values, links):
        if obj_type not in self.store:
            self.store[obj_type] = {}

        self.store[obj_type][unique_id] = {
            "values": dct_values,
            "links": links
        }

    def save(self) :
        raise NotImplemented("This is an abstract class")

    def __getitem__(self, key):
        return self.store[key]

    def __setitem__(self, key, value):
        self.store[key] = value

    def __contains__(self, key):
        return key in self.store
