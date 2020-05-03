import pyGeno.configuration as conf


class ABSObject :
    """

    """
    
    def reset(self):
        self.db = conf.get_database()
        self.collection = self.__class__.__name__

    def _curate(self) :
        pass

    def pygeno_filter_to_aql(self, pygeno_filter, object_name):
        def _rec_rename(pygeno_filter, object_name):
            if pygeno_filter.right_statement:
                _rec_rename(pygeno_filter.statement, object_name)
                _rec_rename(pygeno_filter.right_statement, object_name)
            else :
                pygeno_filter.statement = "%s.%s" %(object_name, pygeno_filter.statement)

        _rec_rename(pygeno_filter, object_name)
        fi = pygeno_filter.to_str().replace("and", "&&").replace("or", "||")
        aql_filters = "FILTER %s" % fi

        aql = """FOR {obj} IN {collection}
        {filters}
        RETURN {obj} 
        """.format(collection=self.collection_name, obj=object_name, filters=aql_filters)

        return aql

    def get(self, type_name, pygeno_filter=None, **lazy_args):
        if type(type_name) is str:
            obj_type = OBJECT_TYPES[type_name]
        else :
            obj_type = type_name

        if pygeno_filter and len(lazy_args) > 0 :
            raise ValueError("You can only use one type of queries at a time")

        if pygeno_filter:
            aql = pygeno_filter_to_aql(pygeno_filter_to_aql, "obj")
            objects = self.db.AQLQuery(aql, batchSize=1000, rawResults=True)
        else :
            objects = self.db[self.collection].fetchByExample(lazy_args, batchSize=1000, rawResults=True)

        for obj in objects;
            yield obj_type(obj)

class Genome(ABSObject) :
    """
    This is the entry point to pyGeno::
        
        myGeno = Genome(name = 'GRCh37.75', SNPs = ['RNA_S1', 'DNA_S1'], SNPFilter = MyFilter)
        for prot in myGeno.get(Protein) :
            print prot.sequence
    
    """

    def _curate(self) :
        self["species"] = self.species.lower()

    def __len__(self) :
        return self["length"]

    def __init__(self, SNPs = None, SNPFilter = None,  *args, **kwargs) :
        pass
        # if type(SNPs) is str :
        #     self.SNPsSets = [SNPs]
        # else :
        #     self.SNPsSets = SNPs

        # # print "pifpasdf", self.SNPsSets

        # if SNPFilter is None :
        #     self.SNPFilter = SF.DefaultSNPFilter()
        # else :
        #     if issubclass(SNPFilter.__class__, SF.SNPFilter) :
        #         self.SNPFilter = SNPFilter
        #     else :
        #         raise ValueError("The value of 'SNPFilter' is not an object deriving from a subclass of SNPFiltering.SNPFilter. Got: '%s'" % SNPFilter)

        # self.SNPTypes = {}
        
        # if SNPs is not None :
        #     f = RabaQuery(SNPMaster, namespace = self._raba_namespace)
        #     for se in self.SNPsSets :
        #         f.addFilter(setName = se, species = self.species)

        #     res = f.run()
        #     if res is None or len(res) < 1 :
        #         raise ValueError("There's no set of SNPs that goes by the name of %s for species %s" % (SNPs, self.species))

        #     for s in res :
        #         # print s.setName, s.SNPType
        #         self.SNPTypes[s.setName] = s.SNPType

    def __str__(self) :
        return "Genome: %s/%s SNPs: %s" %(self.species, self.name, self.SNPTypes)