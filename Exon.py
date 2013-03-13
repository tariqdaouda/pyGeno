import numpy as N
import re, copy, sys, random

from tools import UsefulFunctions as uf
from Protein import Protein

from tools.BinarySequence import NucBinarySequence

class Exon:

	def __init__(self, x1, x2, number, transcript, typ, SNVsFilter = None, startCodon = -1, stopCodon = -1) :
		r"""An exon, the sequence is set according to gene strand, if it's '-' the sequence is the complement.
		A CDS is a couple of coordinates that lies inside of the exon.
		SNVsFilter is a fct that defines wich SNVs are included in the sequence"""
		
		self.type = typ
		self.transcript = transcript
		self.startCodonPos = startCodon
		self.stopCodonPos = stopCodon
		
		xx1, xx2 = int(x1), int(x2)
		if xx1 < xx2 :
			self.x1, self.x2 = xx1, xx2
		else :
			self.x1, self.x2 = xx2, xx1
		
		self.number = int(number)
		self.CDS = None
		
		seq = self.transcript.gene.chromosome.getSequence(x1, x2+1, SNVsFilter)
		if self.transcript.gene.strand == '+' :
			self.sequence = seq
		else :
			self.sequence = uf.reverseComplement(seq)

	def hasCDS(self) :
		if self.CDS != None :
			return True
		return False
	
	def setCDS(self, x1, x2):

		if self.CDS != None and (self.CDS[0] != x1 or self.CDS[1] != x2+1):
			print "==>Warning, Exon.setCDS() : exon %s already has a CDS defined, new CDS: %s " % (self, (x1, x2+1))
		
		xx1, xx2 = int(x1), int(x2)
		if xx1 < xx2 :
			self.CDS = (xx1, xx2)
		else :
			self.CDS = (xx2, xx1)
		
	def getCDSLength(self) :
		return len(self.getCDSSequence())

	
	def getCDSSequence(self) :
		try :
			if self.transcript.gene.strand == '+' :
				return self.sequence[self.CDS[0]-self.x1:self.CDS[1]-self.x1]
			else :
				x1 = self.CDS[0]-self.x1
				x2 = self.CDS[1]-self.x1
				return self.sequence[len(self.sequence)-x2:len(self.sequence)-x1]
		except TypeError:
			return ''
	
	def __getitem__(self, i) :
		return self.sequence[i]
	
	def __str__(self) :
		return """EXON, number: %d, (x1, x2): (%d, %d), cds: %s / %s """ %(self.number, self.x1, self.x2, str(self.CDS), str(self.transcript))
		
	def __len__(self) :
		return len(self.sequence)
