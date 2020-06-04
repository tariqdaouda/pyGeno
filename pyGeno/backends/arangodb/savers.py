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
            else :
                print("TRUNCATING (temporary for tests, should be removed)", colname)
                self.db[colname].truncate()

        # for graph in schemas.ALL_GRAPHS:
        #     graph_name = graph.__name__
        #     if graph_name not in self.db:
        #         self.db.createGraph(graph_name)

    def create_objects(self):
        for colname, objs in self.data.items():
            for key, obj in objs.items():
                try :
                    doc = self.db[colname][key]
                except :
                    doc = self.db[colname].createDocument()
                    doc["_key"] = key
                doc.set(obj)
                doc.save()

    def create_links(self):
        def _get_collection(from_type, to_type):
            if from_type[0] > to_type[0] :
                name = "link_%s_%s" % (from_type, to_type)
            else :
                name = "link_%s_%s" % (to_type, from_type)

            try :
                return self.db[name]
            except :
                return self.db.createCollection(name=name, className="Edges")

        for link in self.links:
            from_key = "%s/%s" % (link["from"]["type"], link["from"]["id"])
            to_key = "%s/%s" % (link["to"]["type"], link["to"]["id"])
            edge = _get_collection(link["from"]["type"], link["to"]["type"]).createEdge()
            edge["_from"] = from_key
            edge["_to"] = to_key
            edge.save()
        
    def save(self) :
        self.init_db()
        self.create_objects()
        self.create_links()