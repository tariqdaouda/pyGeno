class DatabaseConfiguration_ABS:
    """The general interface for objects storing database configurations"""
    def __init__(self, conf=None):
        if conf:
            self.reset(conf)

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
        self.store = {
            "chromosomes": {},
            "genes": {},
            "transcripts": {},
            "proteins": {},
            "exons": {},
        }

    def add(obj_type, unique_id, dct_values):
        self.store[obj_type][unique_id] = dct_values

    def save(self) :
        raise NotImplemented("This is an abstract class")
        # self.conf.db.beginTransaction()
        
        # for c in self.genes.values() :
        #     c.save()
        #     conf.removeFromDBRegistery(c)
            
        # for c in self.transcripts.values() :
        #     c.save()
        #     conf.removeFromDBRegistery(c.exons)
        #     conf.removeFromDBRegistery(c)
        
        # for c in self.proteins.values() :
        #     c.save()
        #     conf.removeFromDBRegistery(c)
        
        # self.conf.db.endTransaction()
        
        # del(self.genes)
        # del(self.transcripts)
        # del(self.proteins)
        # del(self.exons)
        
        # self.genes = {}
        # self.transcripts = {}
        # self.proteins = {}
        # self.exons = {}

        # gc.collect()

    # def save_chros(self) :
    #     pBar = ProgressBar(nbEpochs = len(self.chromosomes))
    #     for c in self.chromosomes.values() :
    #         pBar.update(label = 'Chr %s' % c.number)
    #         c.save()
    #     pBar.close()