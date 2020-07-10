from . import configuration as conf
import pyGeno.tools.UsefulFunctions as uf
# from .pyGenoObjectBases import *

# from .Chromosome import Chromosome
# from .Gene import Gene
# from .Transcript import Transcript
# from .Protein import Protein
# from .Exon import Exon
from . import SNPFiltering as SF
# from .SNP import *

from . backends.rich_query import RichFilter

class GenomeObject:
    """docstring for GenomeObject"""
    def __init__(self, type_name, data=None, unique_id=None):
        self.type_name = type_name
        self.data = data
        self.unique_id = unique_id
        self.fetcher = conf.get_backend().get_query_handler()

    def get(self, object_type, query=None, **lazy_args):
        ret = None
        if (query is not None and not type(query) is dict) and len(lazy_args) > 0:
            raise ValueError("Only dict type queries are compatible with lazy_args")

        if query is None:
            ret = self.fetcher.lazy_query(anchor_type=self.type_name, anchor_id=self.unique_id, fetch_type=object_type)

        if type(query) is dict:
            if len(lazy_args) > 0:
               query = query.update(lazy_args)
            ret = self.fetcher.dict_query(anchor_type=self.type_name, anchor_id=self.unique_id, fetch_type=object_type, dct=query)
        elif len(lazy_args) > 0:
           ret = self.fetcher.lazy_query(anchor_type=self.type_name, anchor_id=self.unique_id, fetch_type=object_type, **lazy_args)

        if type(query) is RichFilter:
           ret = self.fetcher.rich_query(anchor_type=self.type_name, anchor_id=self.unique_id, fetch_type=object_type, pygeno_filter=query)

        if ret is None :
            raise ValueError("Invalid filter type:", query)

        for dct in ret:
            yield GenomeObject(type_name=object_type, data=dct, unique_id=dct["unique_id"])

    @classmethod
    def help(self):
        import json
        fetcher = conf.get_backend().get_query_handler()
        ret = json.dumps(fetcher.get_example("Genome"), indent=4)
        print("Example schema:\n%s" % ret)

    def __getattr__(self, key):
        if key in self.data:
            return self.data[key]

        try:
            res = list(self.get(key.capitalize()))
            if len(res) == 1:
                res = res[0]
        except :
            if key[-1] == "s":
                key = key[:-1]
                res = list(self.get(key.capitalize()))
            else:
                raise KeyError("Unknown attribute: %s" % key)
        return res

    def __getitem__(self, key):
        return self.data[key]

    def __repr__(self) :
        return "%s: %s" %(self.type_name, self.data)

class Transript(object):
    """docstring for Transript"""
    def __init__(self, arg):
        super(Transript, self).__init__()
        self.arg = arg
        
class Genome(GenomeObject) :
    """
    This is the entry point to pyGeno::
        
        myGeno = Genome(name = 'GRCh37.75', SNPs = ['RNA_S1', 'DNA_S1'], SNPFilter = MyFilter)
        for prot in myGeno.get(Protein) :
            print prot.sequence
    
    """
    # _wrapped_class = Genome_Raba

    def __init__(self, name, SNPs = None, SNPFilter = None) :
        super(Genome, self).__init__(type_name="Genome", unique_id=name)

        res_data = list(self.fetcher.lazy_query("Genome", anchor_id=name))
        try:
            self.data = res_data[0]
        except KeyError as e:
            raise KeyError("Can't find Genome: '%s'")
        
        if len(res_data) > 1:
            raise KeyError("Filter is not unique, found %d objects matching" % len(res_data))

        if type(SNPs) is str :
            self.SNPsSets = [SNPs]
        else :
            self.SNPsSets = SNPs

        if SNPFilter is None :
            self.SNPFilter = SF.DefaultSNPFilter()
        else :
            if issubclass(SNPFilter.__class__, SF.SNPFilter) :
                self.SNPFilter = SNPFilter
            else :
                raise ValueError("The value of 'SNPFilter' is not an object deriving from a subclass of SNPFiltering.SNPFilter. Got: '%s'" % SNPFilter)

        self.SNPTypes = {}
        
        if SNPs is not None :
            raise NotImplemented("Not yet implemented")

            f = RabaQuery(SNPMaster, namespace = self._raba_namespace)
            for se in self.SNPsSets :
                f.addFilter(setName = se, species = self.species)

            res = f.run()
            if res is None or len(res) < 1 :
                raise ValueError("There's no set of SNPs that goes by the name of %s for species %s" % (SNPs, self.species))

            for s in res :
                self.SNPTypes[s.setName] = s.SNPType

    def __str__(self) :
        return "Genome: %s/%s SNPs: %s" %(self.species, self.unique_id, self.SNPTypes)
