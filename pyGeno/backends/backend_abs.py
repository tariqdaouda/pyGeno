from . import rich_query

class DatabaseConfiguration_ABS:
    """The general interface for objects storing database configurations"""
    def __init__(self, conf=None):
        if conf:
            self.reset(conf)
        self.saver = None
        self.query_handler = None

    def list_genomes(self):
        raise NotImplemented("This is an abstract class")

    def get_query_handler(self):
        raise NotImplemented("This is an abstract class")

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
        self.data = {}
        self.sequences = {}
        self.links = list()
        self.genome_id = None

    def get_start_end(self, gtf_line):
        """converts from 1 based to 0 based"""
        start = int(gtf_line['start']) - 1
        end = int(gtf_line['end'])
        if start > end :
            start, end = end, start
        return start, end

    def fix_numbers(self, data):
        """make sure int are inte and not strings"""
        if "frame" in data:
            if data["frame"] != ".":
                data["frame"] = int(data["frame"])
        if "exon_number" in data:
            data["exon_number"] = int(data["exon_number"]) 
        return data

    def set_genome(self, unique_id, **kw_data):
        """set the genome"""
        self.data["Genome"] = {"name": unique_id}
        self.data["Genome"][unique_id] = kw_data
        self.genome_id = unique_id

    def _get_link(self, from_type, from_id, to_type, to_id):
        """return a link (edge) connecting two objects"""
        data =  {
            "from": {
                "type": from_type.capitalize(),
                "id": from_id 
            },
            "to":{
                "type": to_type.capitalize(),
                "id": to_id 
            }
        }
        return data

    def link(self, gtf_line, unique_id=None):
        """create a link (edge) between two objects"""
        if unique_id is None:
            unique_id = self.get_uid(gtf_line)
        
        self.links.append(
            self._get_link(self.get_type(gtf_line), unique_id, "Genome", self.genome_id)
        )
        self.links.append(
            self._get_link(self.get_type(gtf_line), unique_id, "Chromosome", gtf_line["seqname"])
        )

        for key in gtf_line["attributes"]:
            if key[-3:] == "_id":
                skey = key.split("_")
                if skey[0] != gtf_line["feature"]:
                    # print(gtf_line["feature"], unique_id, skey[0], gtf_line["attributes"][key])
                    self.links.append(
                        self._get_link(gtf_line["feature"], unique_id, skey[0], gtf_line["attributes"][key])
                    )

    def add(self, gtf_line):
        """add an object to the registery"""
        if self.genome_id is None:
            raise ValueError("Must set genome first")

        if "Chromosome" not in self.data:
            self.data["Chromosome"] = {}

        if gtf_line["seqname"] not in self.data["Chromosome"]: 
            self.data["Chromosome"][gtf_line["seqname"]] = {"number": gtf_line["seqname"]}

        obj_type = self.get_type(gtf_line)
        if obj_type not in self.data:
            self.data[obj_type] = {}

        unique_id = self.get_uid(gtf_line)
        data = gtf_line.to_dct()
        data["start"], data["end"] = self.get_start_end(gtf_line)
        data = self.fix_numbers(data)
        data["genome_id"] = self.genome_id
        self.data[obj_type][unique_id] = data
        self.link(gtf_line)
        self.add_protein(gtf_line)

    def add_protein(self, gtf_line):
        """add a protein object. Proteins are not lines in the gtf, the info has to be extrected from other objects"""
        if "Protein" not in self.data:
            self.data["Protein"] = {}

        if gtf_line["protein_id"] and gtf_line["protein_id"] not in self.data["Protein"]:
            data = {}
            for key, value in gtf_line.to_dct().items():
                if key[-3:] == "_id":
                    data[key] = value
            self.data["Protein"][gtf_line["protein_id"]] = data
            self.link(gtf_line, gtf_line["protein_id"])

    def contains(self, obj_type, unique_id):
        """ask if an bject of of specific types exists in the registery"""
        obj_type = obj_type.capitalize()
        try :
            return self.data[obj_type][unique_id]
        except KeyError:
            return None

    def __contains__(self, key):
        if key in self.data:
            return True
        for typ in self.data:
            if key in self.data[typ]:
                return True
        return False

    def has_sequence(self, key):
        """ask if a sequence is present for a given id"""
        return key.lower() in self.sequences

    def add_sequence(self, key, seq):
        """add a sequence to the registery"""
        self.sequences[key.lower()] = seq

    def get_subsequence(self, seq_id, start, end):
        """get a subsequence from the registery"""
        return self.sequences[seq_id.lower()][start: end]

    def get_type(self, gtf_line):
        """get ovject type name"""
        return gtf_line["feature"].capitalize()

    def get_uid(self, gtf_line):
        """
        makes a unique id based on information in the line. This is useful for things that do not
        already have a unique id, like stop_codons
         """
        if gtf_line["feature"].lower() in ["cds", "exon", "start_codon", "stop_codon"]:
            return "%s_%s" % (gtf_line['transcript_id'], gtf_line['exon_number'])
        elif gtf_line["feature"].lower() in ["utr"]:
            return "utr_" + gtf_line['transcript_id'] + str(gtf_line['start']) + str(gtf_line['end'])

        return gtf_line["attributes"][gtf_line["feature"].lower() + "_id"]

    def save(self) :
        """save registery into backend"""
        raise NotImplemented("This is an abstract class")

    def __getitem__(self, obj_type):
        obj_type = obj_type.capitalize()
        return self.data[obj_type]

class QueryHandler_ABS:
    """
    Saves genome into database
    """
    def __init__(self, database_configuration):
        self.database_configuration = database_configuration

    def get(self, *args, **kwargs):
        if len(args) == 0:
            return self.dict_query(kwargs)
        elif len(args) == 1:
            if type(args[0]) is dict:
                return self.lazy_query(args[0])
            elif type(args[0]) is rich_query.Filter:
                return self.rich_query(args[0])
            else :
                raise ValueError("Unknown query type: %s" % args[0])
        else :
            raise ValueError("Get only takes 1 argument, dict of rich query, or a layze query: a = x, b = y, etc..." % args[0])

    def get_example(self, object_type):
        raise NotImplemented("Must be implemented in child")

    def rich_query(self, pygeno_filter, object_name):
        raise NotImplemented("Must be implemented in child")

    def lazy_query(self, type_name, **lazy_args):
        raise NotImplemented("Must be implemented in child")
        
    def dict_query(self, anchor_type, type_name, dct):
        raise NotImplemented("Must be implemented in child")
