import os, glob, gzip, tarfile, shutil, time, sys, pickle, tempfile
import urllib.request, urllib.error, urllib.parse
from contextlib import closing
from configparser import SafeConfigParser

from pyGeno.tools.ProgressBar import ProgressBar
import pyGeno.configuration as conf

from pyGeno.tools.parsers.GTFTools import GTFFile
from pyGeno.tools.ProgressBar import ProgressBar
from pyGeno.tools.io import printf

def unpack_datawrap(package_file) :
    pFile = tarfile.open(package_file)
    
    package_dir = tempfile.mkdtemp(prefix = "pyGeno_import_")
    if os.path.isdir(package_dir) :
        shutil.rmtree(package_dir)
    os.makedirs(package_dir)

    for mem in pFile :
        pFile.extract(mem, package_dir)

    return package_dir

def get_file(fil, directory) :
    if fil.find("http://") == 0 or fil.find("ftp://") == 0 :
        printf("Downloading file: %s..." % fil)
        final_file = os.path.normpath('%s/%s' %(directory, fil.split('/')[-1]))
        urllib.request.urlretrieve (fil, final_file)
        printf('done.')
    else :
        final_file = os.path.normpath('%s/%s' %(directory, fil))
    
    return final_file

def delete_genome(species, name) :
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

def import_genome(package_file, batch_size = 50, verbose = 0) :
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

    printf('Importing genome package: %s... (This may take a while)' % package_file)

    is_dir = False
    if not os.path.isdir(package_file) :
        package_dir = unpack_datawrap(package_file)
    else :
        is_dir = True
        package_dir = package_file

    parser = SafeConfigParser()
    parser.read(os.path.normpath(package_dir+'/manifest.ini'))
    packageInfos = parser.items('package_infos')

    genome_name = parser.get('genome', 'name')
    species = parser.get('genome', 'species')
    genomeSource = parser.get('genome', 'source')
    
    gtf_file = get_file(parser.get('gene_set', 'gtf'), package_dir)
    chromosome_files = {
        key: get_file(fil, package_dir) for key, fil in parser.items('chromosome_files') 
    }
    
    database = conf.get_backend()
    saver = database.load_saver()

    printf("Importing:\n\t%s\nGenome:\n\t%s\n..."  % (clean_string(packageInfos), clean_string(parser.items('genome'))))
    saver.set_genome(genome_name, species=species, source=genomeSource, packageInfos=packageInfos)
    saver = import_genome_objects(saver, gtf_file, chromosome_files, batch_size, verbose)
    saver.save()
            
    if not is_dir :
        shutil.rmtree(package_dir)
            
def import_genome_objects(saver, gtf_file_path, chromosome_files, batch_size, verbose = 0) :
    """verbose must be an int [0, 4] for various levels of verbosity"""
        
    printf('Importing gene set infos from %s...' % gtf_file_path)
    
    printf("Parsing gene set...")
    gtf = GTFFile(gtf_file_path, gziped = True)
    
    printf('Done. Importation begins!')
    pBar = ProgressBar(nbEpochs = len(gtf))
    for line in gtf :
        chro_fasta_file = chromosome_files[line["seqname"].lower()]
        saver.add(line)
        if not saver.has_sequence(line["seqname"]):
            header, chro_sequence = import_sequence(chro_fasta_file)
            saver["Chromosome"][line["seqname"]]["fasta_header"] = header
            saver.add_sequence(line["seqname"], chro_sequence)

    return saver

def import_sequence(fasta_file) :
    """Serializes fastas file returns original header and sequence as an upper case single line"""
    with gzip.open(fasta_file, 'rt') as f:
        header = f.readline()[:-1]
        sequence = f.read().upper().replace('\n', '').replace('\r', '')
    
        return header, sequence
