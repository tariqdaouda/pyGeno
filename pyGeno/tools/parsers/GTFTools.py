import gzip
import os

class GTFEntry(object) :
    def __init__(self, gtfFile, lineNumber) :
        """A single entry in a GTF file"""
        
        self.lineNumber = lineNumber
        self.gtfFile = gtfFile
        self.data = gtfFile.lines[lineNumber][:-2].split('\t') #-2 remove ';\n'
        proto_atts = self.data[gtfFile.legend['attributes']].strip().split('; ')
        atts = {}
        for a in proto_atts :
            sa = a.split(' ')
            k = sa[0]
            v = sa[1].replace('"', '')
            if k not in atts:
                atts[k] = v
            else:
                if type(atts[k]) != list:
                    atts[k] = [atts[k], v]
                else:
                    atts[k].append(v)
        self.data[gtfFile.legend['attributes']] = atts
    
    def __getitem__(self, k) :
        try :
            return self.data[self.gtfFile.legend[k]]
        except KeyError :
            try :
                return self.data[self.gtfFile.legend['attributes']][k]
            except KeyError :
                #return None
                raise KeyError("Line %d does not have an element %s.\nline:%s" %(self.lineNumber, k, self.gtfFile.lines[self.lineNumber]))
    
    def __repr__(self) :
        return "<GTFEntry line: %d>" % self.lineNumber
    
    def __str__(self) :
        return  "<GTFEntry line: %d, %s>" % (self.lineNumber, str(self.data))

class GTFFile(object) :
    """This is a simple GTF2.2 (Revised Ensembl GTF) parser, see http://mblab.wustl.edu/GTF22.html for more infos"""
    def __init__(self, filename, gziped = False) :
        
        self.filename = filename
        self.legend = {'seqname' : 0, 'source' : 1, 'feature' : 2, 'start' : 3, 'end' : 4, 'score' : 5, 'strand' : 6, 'frame' : 7, 'attributes' : 8}
        self.gziped = gziped

        if gziped : 
            f = gzip.open(filename, 'rt')
        else :
            f = open(filename)
        
        self.lines = []
        for l in f :
            if l[0] != '#' and l != '' :
                self.lines.append(l)
        f.close()
        
        self.currentPos = -1
    
    def get(self, line, elmt) :
        """returns the value of the field 'elmt' of line 'line'"""
        return self[line][elmt]

    def get_transcripts(self, transcript_ids=None):
        """returns genes with its transcripts and associated exons and CDSs from a GTF
        if transcript_ids is used, only these transcripts will be returned"""
        
        if transcript_ids is not None:
            #transcript_ids = [y for x in transcript_ids for y in [x, x.split('.')[0]]]
            transcript_ids = set(transcript_ids) if transcript_ids is not None else transcript_ids
        
        gtf = iter(self)
        
        end_reached = False
        
        # Iterate through all genes in file
        
        # first line
        try:
            line = next(gtf)
        except:
            raise ('Empty GTF file') 
        
        # pop-out lines that come before first gene
        while line['feature'] != 'gene':
            try:
                line = next(gtf)
            except:
                raise ('Never encountered a gene in the GTF') 
        
        # Get all genes
        while line['feature'] == 'gene':
            gene = line
            transcripts = []
            exons = []
            cdss = []
            
            # pop-out all lines that come before first transcript in gene
            while line['feature'] != 'transcript':
                try:
                    line = next(gtf)
                except:
                    raise StopIteration('Last gene %s contains no transcripts' % gene['gene_id']) 
            
            # Get all transcripts in gene
            while line['feature'] == 'transcript':
                tr = line
                ex = []
                cds = []
                
                try:
                    line = next(gtf)
                except:
                    raise StopIteration('Last transcript %s contains no features' % tr['transcript_id']) 
                 
                # Get all features in transcript
                while line['feature'] != 'transcript' and line['feature'] != 'gene':
                    if line['feature'] == 'exon':
                        ex.append(line)
                    if line['feature'] == 'CDS':
                        cds.append(line)
                    
                    try:
                        line = next(gtf)
                    except StopIteration:
                        end_reached = True
                        break
                
                if transcript_ids is None or tr[0]['transcript_id'] in transcript_ids:
                    transcripts.append(tr)
                    exons.append(ex)
                    cdss.append(cds)
                    
            assert len(transcripts)
            yield (gene, transcripts, exons, cdss)
        
        assert end_reached # in a normal GTF file, this value should be true
        assert self.currentPos == len(self)
    
    def gtf2bed(self, bed_filename, feature='transcripts'):
        """Transform gtf to bed6/bed12 and saves the output to file"""
        
        if feature == 'transcripts':
            return self.gtf2bed_transcripts(bed_filename)
        elif feature == 'exons':
            return self.gtf2bed_exons(bed_filename)
        elif feature == 'cds':
            return self.gtf2bed_cds(bed_filename)
        else:
            raise ValueError('Accepted feature types are "transcripts", "exons" and "cds".')
    
    def gtf2bed_transcripts(self, bed_filename):
        """Retrieves transcript information from gtf in bed12 format"""
        
        f_out = open(bed_filename, 'w')
        for gene, transcripts, exons, cdss in self.get_transcripts():
            for tr, ex, cds in zip(transcripts, exons, cdss):
                has_exon = bool(len(ex))
                assert has_exon
                
                chromosome = tr['seqname']
                tid = tr['transcript_id']
                strand = tr['strand']
                start = int(tr['start']) - 1  # 0-base
                end = int(tr['end'])
                
                try:
                    has_stop_codon = 'cds_end_NF' not in set(tr['tag']) and len(cds)
                except KeyError:
                    has_stop_codon = True and len(cds)
                
                exon_lengths = []
                exon_start_pos = []
                
                for ex_line in ex:
                    assert ex_line['transcript_id'] == tid
                    s = int(ex_line['start']) - 1  # 0-base
                    e = int(ex_line['end'])
                    exon_lengths.append(e-s)
                    exon_start_pos.append(s-start)
                    exon_number = ex_line['exon_number']
                
                if len(cds):
                    cds_start_pos = float('inf')
                    cds_end_pos = -1
                    for cds_line in cds:
                        assert cds_line['transcript_id'] == tid
                        s = int(cds_line['start']) - 1  # 0-base
                        e = int(cds_line['end'])
                        cds_start_pos = min(cds_start_pos, s)
                        cds_end_pos = max(cds_end_pos, e)
                else:
                    cds_start_pos = start
                    cds_end_pos = end
                
                if strand == '+':
                    if has_stop_codon:
                        cds_end_pos += 3
                elif strand == '-':
                    if has_stop_codon:
                        cds_start_pos -= 3
                    exon_lengths = exon_lengths[::-1]
                    exon_start_pos = exon_start_pos[::-1]
                else:
                    raise ValueError("Unrecognizable strand value: '%s'." % strand)
                
                entry = [chromosome, start, end, tid, '.', strand, cds_start_pos,
                         cds_end_pos, '.', exon_number,
                         ','.join([str(x) for x in exon_lengths]),
                         ','.join([str(x) for x in exon_start_pos])]
                
                f_out.write('\t'.join([str(x) for x in entry]) + '\n')
        
        f_out.close()
        
        return True
    
    def gtf2bed_exons(self, bed_filename, join_overlaps=True):
        """Retrieves exon information from gtf in bed6 format"""
        
        chromosome_list = []
        exon_positions_dict = dict()
        for gene, transcripts, exons, cdss in self.get_transcripts():
            gene_exon_positions = []
            for tr, ex, cds in zip(transcripts, exons, cdss):
                if len(ex):
                    chromosome = tr['seqname']
                    gid = tr['gene_id']
                    strand = tr['strand']
                   
                    if chromosome not in exon_positions_dict:
                        exon_positions_dict[chromosome] = []
                        chromosome_list.append(chromosome)
                    
                    exon_positions = [(int(e['start'])-1, int(e['end'])) for e in ex]
                    
                    if strand == '+':
                        pass
                    elif strand == '-':
                        exon_positions = exon_positions[::-1]
                    else:
                        raise ValueError("Unrecognizable strand value: '%s'." % strand)
                    
                    gene_exon_positions.extend(exon_positions)
             
            # Join overlapping entries
            if join_overlaps:
                gene_exon_positions = sorted(list(set(gene_exon_positions)))
                gene_exon_positions = self._join_ends(gene_exon_positions)
            
            for start, end in gene_exon_positions:
                exon_positions_dict[chromosome].append((chromosome, start, end, gid, '.', strand))
             
        f_out = open(bed_filename, 'w')
        for chromosome in chromosome_list:
            current_positions = list(set(exon_positions_dict[chromosome]))
            current_positions = sorted(current_positions, key=lambda x: (x[1], x[2], x[3]))
            for entry in current_positions:
                f_out.write('\t'.join([str(x) for x in entry]) + '\n')
        f_out.close()
        
        return True

    def gtf2bed_cds(self, bed_filename, join_overlaps=True):
        """Retrieves CDS information from gtf in bed6 format"""
        
        chromosome_list = []
        cds_positions_dict = dict()
        for gene, transcripts, exons, cdss in self.get_transcripts():
            gene_cds_positions = []
            for tr, ex, cds in zip(transcripts, exons, cdss):
                if len(cds):
                    chromosome = tr['seqname']
                    gid = tr['gene_id']
                    strand = tr['strand']
                    try:
                        has_stop_codon = 'cds_end_NF' not in set(tr['tag'])
                    except KeyError:
                        has_stop_codon = True
                    
                    if chromosome not in cds_positions_dict:
                        cds_positions_dict[chromosome] = []
                        chromosome_list.append(chromosome)
                    
                    cds_positions = [(int(c['start'])-1, int(c['end'])) for c in cds]
                    
                    if strand == '+':
                        if has_stop_codon:
                            cds_positions[-1] = (cds_positions[-1][0], cds_positions[-1][1] + 3)
                    elif strand == '-':
                        if has_stop_codon:
                            cds_positions[-1] = (cds_positions[-1][0] - 3, cds_positions[-1][1])
                        cds_positions = cds_positions[::-1]
                    else:
                        raise ValueError("Unrecognizable strand value: '%s'." % strand)
                    
                    gene_cds_positions.extend(cds_positions)
             
            # Join overlapping entries
            if join_overlaps:
                gene_cds_positions = sorted(list(set(gene_cds_positions)))
                gene_cds_positions = self._join_ends(gene_cds_positions)
            
            for start, end in gene_cds_positions:
                cds_positions_dict[chromosome].append((chromosome, start, end, gid, '.', strand))
             
        f_out = open(bed_filename, 'w')
        for chromosome in chromosome_list:
            current_positions = list(set(cds_positions_dict[chromosome]))
            current_positions = sorted(current_positions, key=lambda x: (x[1], x[2], x[3]))
            for entry in current_positions:
                f_out.write('\t'.join([str(x) for x in entry]) + '\n')
        f_out.close()
        
        return True

    def _join_ends(self, positions, kept=[]):
        if not positions:
            return kept
        elif len(positions) == 1:
            start, end  = positions.pop(0)
            return self._join_ends(positions, kept=kept+[(start, end)])
        else:
            start, end = positions.pop(0)
            start_next, end_next = positions.pop(0)
            while end >= start_next and positions:
                start = min(start, start_next)
                end = max(end, end_next)
                start_next, end_next = positions.pop(0)
            return self._join_ends(positions, kept=kept+[(start, end)])
                
    def __iter__(self) :
        self.currentPos = -1
        return self

    def __next__(self) :
        self.currentPos += 1
        try :
            return GTFEntry(self, self.currentPos)
        except IndexError:
            raise StopIteration

    def __getitem__(self, i) :
        """returns the ith entry"""
        if self.lines[i].__class__ is not GTFEntry :
            self.lines[i] = GTFEntry(self, i)
        return self.lines[i]

    def __repr__(self) :
        return "<GTFFile: %s>" % (os.path.basename(self.filename))

    def __str__(self) :
        return "<GTFFile: %s, gziped: %s, len: %d, currentPosition: %d>" % (os.path.basename(self.filename), self.gziped, len(self), self.currentPos)
        return "<GTFFile: %s, gziped: %s, len: %d>" % (os.path.basename(self.filename), self.gziped, len(self))
    
    def __len__(self) :
        return len(self.lines)
