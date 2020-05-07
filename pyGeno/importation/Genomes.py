import os, glob, gzip, tarfile, shutil, time, sys, pickle, tempfile
import urllib.request, urllib.error, urllib.parse
from contextlib import closing
from configparser import SafeConfigParser

from pyGeno.tools.ProgressBar import ProgressBar
import pyGeno.configuration as conf

# from pyGeno.Genome import *
# from pyGeno.Chromosome import *
# from pyGeno.Gene import *
# from pyGeno.Transcript import *
# from pyGeno.Exon import *
# from pyGeno.Protein import *

from pyGeno.tools.parsers.GTFTools import GTFFile
from pyGeno.tools.ProgressBar import ProgressBar
from pyGeno.tools.io import printf

# import gc
#~ import objgraph

# def backUpDB() :
#     """backup the current database version. automatically called by importGenome(). Returns the filename of the backup"""
#     st = time.ctime().replace(' ', '_')
#     fn = conf.pyGeno_RABA_DBFILE.replace('.db', '_%s-bck.db' % st)
#     shutil.copy2(conf.pyGeno_RABA_DBFILE, fn)

#     return fn

def _decompress_package(packageFile) :
    pFile = tarfile.open(packageFile)
    
    packageDir = tempfile.mkdtemp(prefix = "pyGeno_import_")
    if os.path.isdir(packageDir) :
        shutil.rmtree(packageDir)
    os.makedirs(packageDir)

    for mem in pFile :
        pFile.extract(mem, packageDir)

    return packageDir

def _get_file(fil, directory) :
    if fil.find("http://") == 0 or fil.find("ftp://") == 0 :
        printf("Downloading file: %s..." % fil)
        finalFile = os.path.normpath('%s/%s' %(directory, fil.split('/')[-1]))
        urllib.request.urlretrieve (fil, finalFile)
        printf('done.')
    else :
        finalFile = os.path.normpath('%s/%s' %(directory, fil))
    
    return finalFile

def deleteGenome(species, name) :
    """Removes a genome from the database"""

    printf('deleting genome (%s, %s)...' % (species, name))

    conf.db.beginTransaction()
    objs = []
    allGood = True
    try :
        genome = Genome_Raba(name = name, species = species.lower())
        objs.append(genome)
        pBar = ProgressBar(label = 'preparing')
        for typ in (Chromosome_Raba, Gene_Raba, Transcript_Raba, Exon_Raba, Protein_Raba) :
            pBar.update()
            f = RabaQuery(typ, namespace = genome._raba_namespace)
            f.addFilter({'genome' : genome})
            for e in f.run(gen=True) :
                objs.append(e)
        pBar.close()
        
        pBar = ProgressBar(nbEpochs = len(objs), label = 'deleting objects')
        for e in objs :
            pBar.update()
            e.delete()
        pBar.close()
        
    except KeyError as e :
        #~ printf("\tWARNING, couldn't remove genome form db, maybe it's not there: ", e)
        raise KeyError("\tWARNING, couldn't remove genome form db, maybe it's not there: ", e)
        allGood = False
    printf('\tdeleting folder')
    try :
        shutil.rmtree(conf.getGenomeSequencePath(species, name))
    except OSError as e:
        #~ printf('\tWARNING, Unable to delete folder: ', e)
        OSError('\tWARNING, Unable to delete folder: ', e)
        allGood = False
        
    conf.db.endTransaction()
    return allGood

def importGenome(packageFile, batch_size = 50, verbose = 0) :
    """Import a pyGeno genome package. A genome packages is folder or a tar.gz ball that contains at it's root:

    * gziped fasta files for all chromosomes, or URLs from where them must be downloaded
    
    * gziped GTF gene_set file from Ensembl, or an URL from where it must be downloaded
    
    * a manifest.ini file such as::
    
        [package_infos]
        description = Test package. This package installs only chromosome Y of mus musculus
        maintainer = Tariq Daouda
        maintainer_contact = tariq.daouda [at] umontreal
        version = GRCm38.73

        [genome]
        species = Mus_musculus
        name = GRCm38_test
        source = http://useast.ensembl.org/info/data/ftp/index.html

        [chromosome_files]
        Y = Mus_musculus.GRCm38.73.dna.chromosome.Y.fa.gz / or an url such as ftp://... or http://

        [gene_set]
        gtf = Mus_musculus.GRCm38.73_Y-only.gtf.gz / or an url such as ftp://... or http://

    All files except the manifest can be downloaded from: http://useast.ensembl.org/info/data/ftp/index.html
    
    A rollback is performed if an exception is caught during importation
    
    batch_size sets the number of genes to parse before performing a database save. PCs with little ram like
    small values, while those endowed with more memory may perform faster with higher ones.
    
    Verbose must be an int [0, 4] for various levels of verbosity
    """

    def clean_string(items) :
        s = str(items)
        s = s.replace('[', '').replace(']', '').replace("',", ': ').replace('), ', '\n').replace("'", '').replace('(', '').replace(')', '')
        return s.replace('\n', '\n\t')

    printf('Importing genome package: %s... (This may take a while)' % packageFile)

    isDir = False
    if not os.path.isdir(packageFile) :
        packageDir = _decompress_package(packageFile)
    else :
        isDir = True
        packageDir = packageFile

    parser = SafeConfigParser()
    parser.read(os.path.normpath(packageDir+'/manifest.ini'))
    packageInfos = parser.items('package_infos')

    genome_name = parser.get('genome', 'name')
    species = parser.get('genome', 'species')
    genomeSource = parser.get('genome', 'source')
    
    gtf_file = _get_file(parser.get('gene_set', 'gtf'), packageDir)
    
    chromosome_files = {}
    for key, fil in parser.items('chromosome_files') :
        chromosome_files[key] = _get_file(fil, packageDir)
    
    database = conf.get_backend()
    saver = database.load_saver()

    dct_data = dict(name=genome_name, species=species, source=genomeSource, packageInfos=packageInfos)
    saver.add("Genome", genome_name, dct_data, links = {})
    printf("Importing:\n\t%s\nGenome:\n\t%s\n..."  % (clean_string(packageInfos), clean_string(parser.items('genome'))))

    chros = import_genome_objects(saver, gtf_file, genome_name, batch_size, verbose)
    saver.save()

    if False:
        # os.makedirs(seqTargetDir)
        startChro = 0
        pBar = ProgressBar(nbEpochs = len(chros))
        for chro in chros :
            pBar.update(label = "Importing DNA, chro %s" % chro.number)
            length = _importSequence(chro, chromosome_files[chro.number.lower()], seqTargetDir)
            chro.start = startChro
            chro.end = startChro+length
            startChro = chro.end
            chro.save()
        pBar.close()
        
        if not isDir :
            shutil.rmtree(packageDir)
    
    return True

def get_start_end(line):
    start = int(line['start']) - 1
    end = int(line['end'])
    if start > end :
        start, end = end, start
    return start, end
        
def parse_chromosome(saver, line, genome_id):
    number = line['seqname']
    saver.add("Chromosome", number, {"number": number}, links = {"Genome": [genome_id]})

def parse_gene(saver, line, genome_id):
    if line['gene_id'] is not None :
        start, end = get_start_end(line)

        if ("Gene" not in saver) or (line["gene_id"] not in saver["Gene"]):
            dct_data = dict(
                id=line['gene_id'],
                name=line['gene_name'],
                strand=line['strand'],
                biotype=line['gene_biotype'],
                start=start,
                end=end
            )
            links = {
                "Genome": [genome_id],
                "Chromosome": [line['seqname']]
            }
            saver.add("Gene", line['gene_id'], dct_data, links)
        else: 
            if start < saver["Gene"][line['gene_id']]["values"]["start"]:
                saver["Gene"][line['gene_id']]["values"]["start"] = start
            if end > saver["Gene"][line['gene_id']]["values"]["end"]:
                saver["Gene"][line['gene_id']]["values"]["end"] = end

def parse_transcript(saver, line, genome_id):
    try :
        trans_id = line['transcript_id']
        trans_name = line['transcript_name']
        try :
            transcript_biotype = line['transcript_biotype']
        except KeyError :
            transcript_biotype = None
        saver.add(
            "Transcript",
            trans_id,
            {
                "name": line['transcript_name'],
                "biotype": line['transcript_biotype']
            }
        )
    except KeyError :
        trans_id = None
        trans_name = None
        if verbose > 2 :
            printf('\t\tWarning: no transcript_id, name found in line %s' % gtf[i])
            
def import_genome_objects(saver, gtf_file_path, genome, batch_size, verbose = 0) :
    """verbose must be an int [0, 4] for various levels of verbosity"""
        
    printf('Importing gene set infos from %s...' % gtf_file_path)
    
    printf("Parsing gene set...")
    gtf = GTFFile(gtf_file_path, gziped = True)
    
    printf('Done. Importation begins!')

    pBar = ProgressBar(nbEpochs = len(gtf))
    for line in gtf :
        # pBar.update()
        
        # strand = line['strand']
        # gene_biotype = line['gene_biotype']
        # regionType = line['feature']
        # frame = line['frame']
 
        # chro_number = chroN.upper()
        parse_chromosome(saver, line, genome)
        parse_gene(saver, line, genome)
        parse_gene(saver, line, genome)
        
        if False:
            try :
                transId = line['transcript_id']
                transName = line['transcript_name']
                try :
                    transcript_biotype = line['transcript_biotype']
                except KeyError :
                    transcript_biotype = None
            except KeyError :
                transId = None
                transName = None
                if verbose > 2 :
                    printf('\t\tWarning: no transcript_id, name found in line %s' % gtf[i])
            
            if transId is not None :
                if transId not in store.transcripts :
                    if verbose > 1 :
                        printf('\t\tTranscript %s, %s...' % (transId, transName))
                    store.transcripts[transId] = Transcript_Raba()
                    store.transcripts[transId].set(genome = genome, id = transId, chromosome = store.chromosomes[chro_number], gene = store.genes.get(geneId, None), name = transName, biotype=transcript_biotype)
                if store.transcripts[transId].start is None or start < store.transcripts[transId].start:
                    store.transcripts[transId].start = start
                if store.transcripts[transId].end is None or end > store.transcripts[transId].end:
                    store.transcripts[transId].end = end
            
                try :
                    protId = line['protein_id']
                except KeyError :
                    protId = None
                    if verbose > 2 :
                        printf('Warning: no protein_id found in line %s' % gtf[i])

                # Store selenocysteine positions in transcript
                if regionType == 'Selenocysteine':
                    store.transcripts[transId].selenocysteine.append(start)
                        
                if protId is not None and protId not in store.proteins :
                    if verbose > 1 :
                        printf('\t\tProtein %s...' % (protId))
                    store.proteins[protId] = Protein_Raba()
                    store.proteins[protId].set(genome = genome, id = protId, chromosome = store.chromosomes[chro_number], gene = store.genes.get(geneId, None), transcript = store.transcripts.get(transId, None), name = transName)
                    store.transcripts[transId].protein = store.proteins[protId]

                try :
                    exonNumber = int(line['exon_number']) - 1
                    exonKey = (transId, exonNumber)
                except KeyError :
                    exonNumber = None
                    exonKey = None
                    if verbose > 2 :
                        printf('Warning: no exon number or id found in line %s' % gtf[i])
                
                if exonKey is not None :
                    if verbose > 3 :
                        printf('\t\t\texon %s...' % (exonId))
                    
                    if exonKey not in store.exons and regionType == 'exon' :
                        store.exons[exonKey] = Exon_Raba()
                        store.exons[exonKey].set(genome = genome, chromosome = store.chromosomes[chro_number], gene = store.genes.get(geneId, None), transcript = store.transcripts.get(transId, None), protein = store.proteins.get(protId, None), strand = strand, number = exonNumber, start = start, end = end)
                        store.transcripts[transId].exons.append(store.exons[exonKey])
                    
                    try :
                        store.exons[exonKey].id = line['exon_id']
                    except KeyError :
                        pass
                    
                    if regionType == 'exon' :
                        if store.exons[exonKey].start is None or start < store.exons[exonKey].start:
                            store.exons[exonKey].start = start
                        if store.exons[exonKey].end is None or end > store.transcripts[transId].end:
                            store.exons[exonKey].end = end
                    elif regionType == 'CDS' :
                        store.exons[exonKey].CDS_start = start
                        store.exons[exonKey].CDS_end = end
                        store.exons[exonKey].frame = frame
                    elif regionType == 'stop_codon' :
                        if strand == '+' :
                            if store.exons[exonKey].CDS_end != None :
                                store.exons[exonKey].CDS_end += 3
                                if store.exons[exonKey].end < store.exons[exonKey].CDS_end :
                                    store.exons[exonKey].end = store.exons[exonKey].CDS_end
                                if store.transcripts[transId].end < store.exons[exonKey].CDS_end :
                                    store.transcripts[transId].end = store.exons[exonKey].CDS_end
                                if store.genes[geneId].end < store.exons[exonKey].CDS_end :
                                    store.genes[geneId].end = store.exons[exonKey].CDS_end
                        if strand == '-' :
                            if store.exons[exonKey].CDS_start != None :
                                store.exons[exonKey].CDS_start -= 3
                                if store.exons[exonKey].start > store.exons[exonKey].CDS_start :
                                    store.exons[exonKey].start = store.exons[exonKey].CDS_start
                                if store.transcripts[transId].start > store.exons[exonKey].CDS_start :
                                    store.transcripts[transId].start = store.exons[exonKey].CDS_start
                                if store.genes[geneId].start > store.exons[exonKey].CDS_start :
                                    store.genes[geneId].start = store.exons[exonKey].CDS_start
            pBar.close()
            
            store.batch_save()
            
            conf.db.beginTransaction()
            printf('almost done saving chromosomes...')
            store.save_chros()
            
            printf('saving genome object...')
            genome.save()
            conf.db.endTransaction()
            
            conf.db.beginTransaction()
            printf('restoring core indexes...')
            # Genome.ensureGlobalIndex(('name', 'species'))
            # Chromosome.ensureGlobalIndex('genome')
            # Gene.ensureGlobalIndex('genome')
            # Transcript.ensureGlobalIndex('genome')
            # Protein.ensureGlobalIndex('genome')
            # Exon.ensureGlobalIndex('genome')
            Transcript.ensureGlobalIndex('exons')
            
            printf('commiting changes...')
            conf.db.endTransaction()
            
            conf.db.beginTransaction()
            printf('restoring user indexes')
            pBar = ProgressBar(label = "restoring", nbEpochs = len(indexes))
            for idx in indexes :
                pBar.update()
                conf.db.execute(idx[-1].replace('CREATE INDEX', 'CREATE INDEX IF NOT EXISTS'))
            pBar.close()
            
            printf('commiting changes...')
            conf.db.endTransaction()
            
            return list(store.chromosomes.values())

#~ @profile
def _importSequence(chromosome, fastaFile, targetDir) :
    "Serializes fastas into .dat files"

    f = gzip.open(fastaFile, 'rt')
    header = f.readline()
    strRes = f.read().upper().replace('\n', '').replace('\r', '')
    f.close()

    fn = '%s/chromosome%s.dat' % (targetDir, chromosome.number)
    f = open(fn, 'w')
    f.write(strRes)
    f.close()
    chromosome.dataFile = fn
    chromosome.header = header
    return len(strRes)
