import pyArango.collection as COL
import pyArango.graph as GR
import pyArango.validation as VAL

class GenomicLink(COL.Edges):
    _fields = {}

# class GenomeGraph(GR.Graph) :
#     _edgeDefinitions = [
#       GR.EdgeDefinition(
#         "GenomicLink",
#         fromCollections=[
#           "Genome", "Chromosome", "Gene", "Exon", "Transcript", "Protein", "Cds", "Start_codon", "Stop_codon", "Utr"
#         ],
#         toCollections=[
#           "Genome", "Chromosome", "Gene", "Exon", "Transcript", "Protein", "Cds", "Start_codon", "Stop_codon", "Utr"
#         ]
#       )
#     ]
#     _orphanedCollections = []

class Genome(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

    # _fields = {
    #     'name': COL.Field(validators=[]),
    #     'species': COL.Field(validators=[]),
    #     'source': COL.Field(validators=[]),
    #     'packageInfos': COL.Field(validators=[]),
    #     'length': COL.Field()
    # }

class Chromosome(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }
  
    # _fields = {
    #     # 'header': COL.Field(validators=[]),
    #     'number': COL.Field(validators=[]),
    #     'start': COL.Field(),
    #     'end': COL.Field(),
    #     'length': COL.Field()
    # }

class Gene(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }
  
    # _fields = {
    #     "id": COL.Field(validators=[]),
    #     "name": COL.Field(validators=[]),
    #     "strand": COL.Field(validators=[]),
    #     "biotype": COL.Field(validators=[]),    
    #     "start": COL.Field(validators=[]),
    #     "end": COL.Field(validators=[]),
    # }

class Exon(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }
  
    # _fields = {
    #     "id": COL.Field(validators=[]),
    #     "number": COL.Field(validators=[]),
    #     "start": COL.Field(validators=[]),
    #     "end": COL.Field(validators=[]),
    #     "length": COL.Field(validators=[]),
    #     "CDS_length": COL.Field(validators=[]),
    #     "CDS_start": COL.Field(validators=[]),
    #     "CDS_end": COL.Field(validators=[]),
    #     "frame": COL.Field(validators=[]),
    #     "strand": COL.Field(validators=[]),
    # }

class Transcript(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }
  
    # _fields = {
    #     "id": COL.Field(validators=[]),
    #     "name": COL.Field(validators=[]),
    #     "length": COL.Field(),
    #     "start": COL.Field(validators=[]),
    #     "end": COL.Field(validators=[]),
    #     "coding": COL.Field(),
    #     "biotype": COL.Field(),
    #     "selenocysteine": COL.Field(),
    # }

class Protein(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }
  
    # _fields = {
    #     "id": COL.Field(validators=[]),
    #     "name": COL.Field(validators=[]),
    # }

class Cds(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

class Start_codon(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

class Stop_codon(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

class Utr(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

ALL_COLLECTIONS = [GenomicLink, Genome, Chromosome, Gene, Exon, Transcript, Protein, Cds, Start_codon, Stop_codon, Utr]
# ALL_GRAPHS = [GenomeGraph]
