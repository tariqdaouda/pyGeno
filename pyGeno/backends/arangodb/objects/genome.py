import pyGeno.configuration as conf

import pyArango.collection as COL
import pyArango.validation as VAL

from . import objects_abs

class Genome(objects_abs.ABSObject) :
    """
    This is the entry point to pyGeno::
        
        myGeno = Genome(name = 'GRCh37.75', SNPs = ['RNA_S1', 'DNA_S1'], SNPFilter = MyFilter)
        for prot in myGeno.get(Protein) :
            print prot.sequence
    
    """
    
    _properties = {
      "keyOptions" : {
          "allowUserKeys": False,
          "type": "autoincrement",
          "increment": 1,
          "offset": 0,
      }
    }

    _validation = {
        'on_save': True,
        'on_set': True,
        'allow_foreign_fields': True
    }

    _fields = {
        'name': COL.Field(validators=[VAL.NotNull()]),
        'species': COL.Field(validators=[VAL.NotNull()]),
        'source': COL.Field(validators=[VAL.NotNull()]),
        'packageInfos': COL.Field(validators=[VAL.NotNull()]),
        'length': COL.Field()
    }

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
