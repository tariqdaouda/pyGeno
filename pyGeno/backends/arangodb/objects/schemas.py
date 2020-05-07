import pyArango.collection as COL
import pyArango.graph as GR
import pyArango.validation as VAL

class GenomicLink(COL.Edges):
    _fields = {}

class GenomeGraph(GR.Graph) :
    _edgeDefinitions = [
      GR.EdgeDefinition(
        "GenomicLink",
        fromCollections=[
          "Genome", "Chromosome", "Gene", "Exon"
        ],
        toCollections=[
          "Genome", "Chromosome", "Gene", "Exon"
        ]
      )
    ]
    _orphanedCollections = []

class Genome(COL.Collection) :
    _validation = {
        'on_save': True,
        'on_set': True,
        'allow_foreign_fields': False
    }

    _fields = {
        'name': COL.Field(validators=[VAL.NotNull()]),
        'species': COL.Field(validators=[VAL.NotNull()]),
        'source': COL.Field(validators=[VAL.NotNull()]),
        'packageInfos': COL.Field(validators=[VAL.NotNull()]),
        'length': COL.Field()
    }

class Chromosome(COL.Collection) :
    _validation = {
        'on_save': True,
        'on_set': True,
        'allow_foreign_fields': False
    }
  
    _fields = {
        # 'header': COL.Field(validators=[VAL.NotNull()]),
        'number': COL.Field(validators=[VAL.NotNull()]),
        'start': COL.Field(),
        'end': COL.Field(),
        'length': COL.Field()
    }

class Gene(COL.Collection) :
    _validation = {
        'on_save': True,
        'on_set': True,
        'allow_foreign_fields': False
    }
  
    _fields = {
        "id": COL.Field(validators=[VAL.NotNull()]),
        "name": COL.Field(validators=[VAL.NotNull()]),
        "strand": COL.Field(validators=[VAL.NotNull()]),
        "biotype": COL.Field(validators=[VAL.NotNull()]),    
        "start": COL.Field(validators=[VAL.NotNull()]),
        "end": COL.Field(validators=[VAL.NotNull()]),
    }

class Exon(COL.Collection) :
    _validation = {
        'on_save': True,
        'on_set': True,
        'allow_foreign_fields': False
    }
  
    _fields = {
        "id": COL.Field(validators=[VAL.NotNull()]),
        "number": COL.Field(validators=[VAL.NotNull()]),
        "start": COL.Field(validators=[VAL.NotNull()]),
        "end": COL.Field(validators=[VAL.NotNull()]),
        "length": COL.Field(validators=[VAL.NotNull()]),
        "CDS_length": COL.Field(validators=[VAL.NotNull()]),
        "CDS_start": COL.Field(validators=[VAL.NotNull()]),
        "CDS_end": COL.Field(validators=[VAL.NotNull()]),
        "frame": COL.Field(validators=[VAL.NotNull()]),
        "strand": COL.Field(validators=[VAL.NotNull()]),
    }

class Transcript(COL.Collection) :
    _validation = {
        'on_save': True,
        'on_set': True,
        'allow_foreign_fields': False
    }
  
    _fields = {
        "id": COL.Field(validators=[VAL.NotNull()]),
        "name": COL.Field(validators=[VAL.NotNull()]),
        "length": COL.Field(validators=[VAL.NotNull()]),
        "start": COL.Field(validators=[VAL.NotNull()]),
        "end": COL.Field(validators=[VAL.NotNull()]),
        "coding": COL.Field(validators=[VAL.NotNull()]),
        "biotype": COL.Field(validators=[VAL.NotNull()]),
        "selenocysteine": COL.Field(),
    }

class Protein(COL.Collection) :
    _validation = {
        'on_save': True,
        'on_set': True,
        'allow_foreign_fields': False
    }
  
    _fields = {
        "id": COL.Field(validators=[VAL.NotNull()]),
        "name": COL.Field(validators=[VAL.NotNull()]),
    }


ALL_COLLECTIONS = [GenomicLink, Genome, Chromosome, Gene, Exon, Transcript, Protein]
ALL_GRAPHS = [GenomeGraph]
