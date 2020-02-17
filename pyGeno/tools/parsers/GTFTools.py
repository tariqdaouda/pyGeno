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

	def get_transcripts(self, transcript_ids=None):
		"""returns transcript and their associated exons and CDS from GTF (all if transcript_ids is None)"""
		
		if transcript_ids is not None:
			#transcript_ids = [y for x in transcript_ids for y in [x, x.split('.')[0]]]
			transcript_ids = set(transcript_ids) if transcript_ids is not None else transcript_ids
		
		gtf = iter(self)
		line = next(gtf)
		while self.currentPos != len(self):
			try:
				tr = []
				ex = []
				cds = []
				if line['feature'] == 'transcript':
					tr = [line]
					line = next(gtf)
					while line['feature'] != 'transcript':
						if line['feature'] == 'exon':
							ex.append(line)
						if line['feature'] == 'CDS':
							cds.append(line)
						line = next(gtf)
					if transcript_ids is None or tr[0]['transcript_id'] in transcript_ids:
					    yield (tr, ex, cds)
				else:
					line = next(gtf)
			except StopIteration:
				assert tr  # in a normal GTF file, this value should never be empty
				if transcript_ids is None or tr[0]['transcript_id'] in transcript_ids:
					yield (tr, ex, cds)
	
	def gtf2bed(self, bed_filename):
		"""Transform gtf to bed12 abd saves to file"""
		
		f_out = open(bed_filename, 'w')
		for tr, ex, cds in self.get_transcripts():
			tr = tr[0]
			chromosome = tr['seqname']
			tid = tr['transcript_id']
			strand = tr['strand']
			start = int(tr['start']) - 1  # 0-base
			end = int(tr['end'])
			
			has_exon = bool(len(ex))
			try:
				has_stop_codon = 'cds_end_NF' not in set(tr['tag'])
			except KeyError:
				has_stop_codon = True
			
			exon_lengths = []
			exon_start_pos = []
			cds_start_pos = None
			cds_end_pos = None
			
			assert has_exon
			
			for ex_line in ex:
				assert ex_line['transcript_id'] == tid
				s = int(ex_line['start']) - 1  # 0-base
				e = int(ex_line['end'])
				exon_lengths.append(e-s)
				exon_start_pos.append(s-start)
				exon_number = ex_line['exon_number']
			
			for cds_line in cds:
				assert cds_line['transcript_id'] == tid
				s = int(cds_line['start']) - 1  # 0-base
				e = int(cds_line['end'])
				if cds_start_pos is None:
					cds_start_pos = s
				else:
					cds_start_pos = min(cds_start_pos, s)
				if cds_end_pos is None:
					cds_end_pos = e
				else:
					cds_end_pos = max(cds_end_pos, e)
			
			if strand == '+':
				if  cds_start_pos is None:
					cds_start_pos = start
				
				if  cds_end_pos is None:
					cds_end_pos = end
				elif has_stop_codon:
					cds_end_pos = cds_end_pos + 3
			
			elif strand == '-':
				exon_lengths = exon_lengths[::-1]
				exon_start_pos = exon_start_pos[::-1]
				
				if  cds_start_pos is None:
					cds_start_pos = start
				elif has_stop_codon:
				    cds_start_pos = cds_start_pos - 3
				
				if  cds_end_pos is None:
					cds_end_pos = end
			else:
				raise ValueError("Unrecognizable strand value: '%s'." % strand)
			
			entry = [chromosome, start, end, tid, '.', strand, cds_start_pos,
					 cds_end_pos, '.', exon_number,
					 ','.join([str(x) for x in exon_lengths]),
					 ','.join([str(x) for x in exon_start_pos])]
			
			f_out.write('\t'.join([str(x) for x in entry]) + '\n')
		
		f_out.close()
		return True
	
	def get(self, line, elmt) :
		"""returns the value of the field 'elmt' of line 'line'"""
		return self[line][elmt]

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
