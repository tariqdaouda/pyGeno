from pyGenoObject import *

#from rabaDB.setup import *
#RabaConfiguration(conf.pyGeno_RABA_NAMESPACE, conf.pyGeno_RABA_DBFILE)
#from rabaDB.Raba import *
#from rabaDB.filters import RabaQuery
import rabaDB.fields as rf

#from tools import UsefulFunctions as uf
#from Protein import Protein

from tools.BinarySequence import NucBinarySequence

class Exon(pyGenoObject):
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	
	id = rf.Primitive()
	number = rf.Primitive()
	x1 = rf.Primitive()
	x2 = rf.Primitive()
	length = rf.Primitive()
	CDS_x1 = rf.Primitive()
	CDS_x2 = rf.Primitive()
	
	genome = rf.RabaObject('Genome')
	chromosome = rf.RabaObject('Chromosome')
	gene = rf.RabaObject('Gene')
	transcript = rf.RabaObject('Transcript')
	strand = rf.Primitive()
	
	_raba_uniques = [('genome', 'id')]
	
	def __init__(self, importing = False) :
		r"""An exon, the sequence is set according to gene strand, if it's '-' the sequence is the complement.
		A CDS is a couple of coordinates that lies inside of the exon.
		SNVsFilter is a fct that defines wich SNVs are included in the sequence.
		I expect [x1, x2[ (python) not something like [x1, x2](ensembl format)"""
		if importing :
			self.sequence = ''
		else :
			if self.x1 != None and self.x2 != None :
				#xx1, xx2 = int(self.x1), int(self.x2)+1
				xx1, xx2 = int(self.x1), int(self.x2)
				if xx1 < xx2 :
					self.x1, self.x2 = xx1, xx2
				else :
					self.x1, self.x2 = xx2, xx1
			
			if self.number != None :
				self.number = int(self.number)
			
			seq = self.transcript.gene.chromosome.getSequence(x1, x2+1, SNVsFilter)
			if self.transcript.gene.strand == '+' :
				self.sequence = seq
			else :
				self.sequence = uf.reverseComplement(seq)
			
			self.cdsSequence = ''
		
	def save(self) :
		if  self.x2 != None and self.x1 != None :
			self.length = self.x2-self.x1
		if self.number != None :
			self.number = int(self.number)
		
		pyGenoObject.save(self)
	
	def hasCDS(self) :
		if self.CDS_x1 != None and self.CDS_x2 != None:
			return True
		return False
	
	def _getCDSSequence(self) :
		try :
			if self.transcript.gene.strand == '+' :
				return self.sequence[self.CDS_x1-self.x1:self.CDS_x2-self.x1]
			else :
				x1 = self.CDS_x1-self.x1
				x2 = self.CDS_x2-self.x1
				return self.sequence[len(self.sequence)-x2:len(self.sequence)-x1]
		except TypeError:
			return ''
	
	def setCDS(self, x1, x2):
		"Beware! i expect [x1, x2[ (python) not something like [x1, x2](ensembl format)"
		
		#xx1, xx2 = int(x1), int(x2)+1
		xx1, xx2 = int(x1), int(x2)
		if self.CDS_x1 != None and self.CDS_x2 != None and (self.CDS_x1 != xx1 or self.CDS_x2 != xx2):
			print "==>Warning, Exon.setCDS() : exon %s already has a CDS defined, new CDS: %s " % (self, (xxx1, xx2))
		
		if xx1 < xx2 :
			self.CDS_x1, self.CDS_x2 = xx1, xx2
		else :
			self.CDS_x1, self.CDS_x2 = xx2, xx1
	
		self.cdsSequence = self._getCDSSequence()
		
	def getCDSLength(self) :
		return len(self.getCDSSequence())
	
	def pluck(self) :
		"""Returns a plucked object. Plucks the exon off the tree, set the value of self.transcript into str(self.transcript). This effectively disconnects the object and
		makes it much more lighter in case you'd like to pickle it"""
		e = copy.copy(self)
		e.transcript = str(self.transcript)
		return e
		
	def __getitem__(self, i) :
		return self.sequence[i]
	
	def __str__(self) :
		return """EXON, number: %s, (x1, x2): (%s, %s), cds: (%s, %s) > %s """ %(self.number, self.x1, self.x2, self.CDS_x1, self.CDS_x2, str(self.transcript))
		
	def __len__(self) :
		return len(self.sequence)
