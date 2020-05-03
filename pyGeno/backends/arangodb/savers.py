from ..backend_abs import GenomeSaver_ABS
from pyGeno.configuration import system_message
from .objects import schemas

class GenomeSaver(GenomeSaver_ABS):
    """
    Saves genome into database
    """
    def __init__(self, database_configuration):
        super(GenomeSaver, self).__init__(database_configuration)
        self.db = self.database_configuration.database

    def init_db(self):
        for col in schemas.ALL_COLLECTIONS:
            colname = col.__name__
            if colname not in self.db:
                self.db.createCollection(colname)

        for graph in schemas.ALL_GRAPHS:
            graph_name = graph.__name__
            if graph_name not in self.db:
                self.db.createGraph(graph_name)

    def create_objects(self):
        for colname, objs in self.store.items():
            for key, obj in objs.items():
                try :
                    doc = self.db[colname][key]
                except :
                    doc = self.db[colname].createDocument()
                    doc["_key"] = key
                doc.set(obj["values"])
                doc.save()

    def create_links(self):
        for colname, objs in self.store.items():
            for key, obj in objs.items():
                obj_key = "%s/%s" % (colname, key)
                for linked_col, links in obj["links"].items():
                    for link in links :
                        linked_key = "%s/%s" % (linked_col, link)
                        self.db.graphs["GenomeGraph"].link('GenomicLink', obj_key, linked_key, {})

    def save(self) :
        self.init_db()
        self.create_objects()
        self.create_links()