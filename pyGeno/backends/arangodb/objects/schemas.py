import pyArango.collection as COL
import pyArango.graph as GR
import pyArango.validation as VAL

class Genome(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

class Chromosome(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

class Gene(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

class Exon(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

class Transcript(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

class Protein(COL.Collection) :
    _validation = {
        'on_save': False,
        'on_set': False,
        'allow_foreign_fields': True
    }

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

ALL_COLLECTIONS = [Genome, Chromosome, Gene, Exon, Transcript, Protein, Cds, Start_codon, Stop_codon, Utr]
