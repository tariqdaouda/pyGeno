from ..backend_abs import GenomeSaver_ABS
from pyGeno.configuration import system_message
from .objects import schemas
import pyArango
class GenomeSaver(GenomeSaver_ABS):
    """
    Saves genome into database
    """
    def __init__(self, database_configuration):
        super(GenomeSaver, self).__init__(database_configuration)
        self.db = self.database_configuration.database
        self.accepted_contigs = ["Cds", "Ccds", "Exon", "Start_codon", "Stop_codon", "Utr", "Gene", "Transcript"]

    def init_db(self):
        for col in schemas.ALL_COLLECTIONS:
            colname = col.__name__
            if not self.db.hasCollection(str(colname)):
                self.db.createCollection(colname)
            else :
                print("TRUNCATING (temporary for tests, should be removed)", colname)
                self.db[colname].truncate()
    
    def create_objects(self):
        for colname, objs in self.data.items():
            i = 0
            c = 0
            docs = []
            doc = {}
            for key, obj in objs.items():
                doc = self.db[colname].createDocument_(obj)
                doc['_key'] = key
                if colname in self.accepted_contigs:
                    doc["contig"] = self.get_subsequence(doc["seqname"], doc["start"], doc["end"])
                docs.append(doc)
                i = i + 1
                c = c + 1
                if i == len(objs) or i == 1000 or c + i == len(objs):
                    self.db[colname].bulkSave(docs, onDuplicate='ignore')
                    docs = []
                    i = 0

    def create_links(self):
        def _get_collection(from_type, to_type):
            if from_type[0] > to_type[0] :
                name = "%s_X_%s" % (from_type, to_type)
            else :
                name = "%s_X_%s" % (to_type, from_type)

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
        #self.create_links()
