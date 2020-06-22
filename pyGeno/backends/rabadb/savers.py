import rabaDB.rabaSetup
import rabaDB.Raba
from ..backend_abs import GenomeSaver_ABS
from pyGeno.database_configurationiguration import system_message

import gc

class GenomeSaver(GenomeSaver_ABS):
    """
    Saves genome into database
    """
    def __init__(self, database_configuration):
        super(GenomeSaver, self).__init__(database_configuration)
 
    def add(obj_type, unique_id, dct_values):
        self.store[obj_type][unique_id] = dct_values

    def _save_genome(self):
        genome = Genome_Raba()
        genome.set(name = genomeName, species = species, source = genomeSource, packageInfos = packageInfos)
 
    def _save_objects(self):
        self.database_configuration.db.beginTransaction()
        system_message("Saving objects")
        for gene in self.genes.values() :
            gene.save()
            conf.removeFromDBRegistery(gene)
            
        for trans in self.transcripts.values() :
            trans.save()
            conf.removeFromDBRegistery(trans.exons)
            conf.removeFromDBRegistery(trans)
        
        for prot in self.proteins.values() :
            prot.save()
            conf.removeFromDBRegistery(prot)
        
        system_message('commiting changes...')
        self.database_configuration.db.endTransaction()        
        self.store = None
        gc.collect()

    def _save_chormosomes(self):
        pass

    def _drop_indexes(self):
        indexes = self.database_configuration.db.getIndexes()
        
        Genome_Raba.flushIndexes()
        Chromosome_Raba.flushIndexes()
        Gene_Raba.flushIndexes()
        Transcript_Raba.flushIndexes()
        Protein_Raba.flushIndexes()
        Exon_Raba.flushIndexes()
        
        return indexes

    def _restore_indexes(self, indexes):
        self.database_configuration.db.beginTransaction()
        system_message('restoring core indexes...')
        Transcript.ensureGlobalIndex('exons')
        system_message('commiting changes...')
        self.database_configuration.db.endTransaction()

        self.database_configuration.db.beginTransaction()
        system_message('restoring user indexes')
        pBar = ProgressBar(label = "restoring", nbEpochs = len(indexes))
        for idx in indexes :
            pBar.update()
            self.database_configuration.db.execute(idx[-1].replace('CREATE INDEX', 'CREATE INDEX IF NOT EXISTS'))
        pBar.close()
        
        system_message('commiting changes...')
        self.database_configuration.db.endTransaction()

    def save(self) :
        system_message('Backuping indexes...')

        system_message("Dropping all your indexes (don't worry I'll restore them later)...")
        indexes = self._drop_indexes()
        self._save_objects()
        self._restore_indexes(indexes)
