from pyGeno.Genome import Genome
from pyGeno.Gene import Gene
from pyGeno.Transcript import Transcript
from pyGeno.Protein import Protein
from pyGeno.Exon import Exon

from pyGeno.SNPFiltering import SNPFilter
from pyGeno.SNPFiltering import SequenceSNP

def printing(gene) :
	print "printing reference sequences\n-------"
	for trans in gene.get(Transcript) :
		print "\t-Transcript name:", trans.name
		print "\t-Protein:", trans.protein.sequence
		print "\t-Exons:"
		for e in trans.exons :
			print "\t\t Number:", e.number
			print "\t\t-CDS:", e.CDS
			print "\t\t-Strand:", e.strand
			print "\t\t-CDS_start:", e.CDS_start
			print "\t\t-CDS_end:", e.CDS_end

def printVs(refGene, presGene) :
	print "Vs personalized sequences\n------"

	#iterGet returns an iterator instead of a list (faster)
	for trans in presGene.iterGet(Transcript) :
		refProt = refGene.get(Protein, id = trans.protein.id)[0]
		persProt = trans.protein
		print persProt.id
		print "\tref:" + refProt.sequence[:20] + "..."
		print "\tper:" + persProt.sequence[:20] + "..."
		print

def fancyExonQuery(gene) :
	e = gene.get(Exon, {'CDS_start >' : 2655029, 'CDS_end <' : 2655200})[0]
	print "An exon with a CDS in ']2655029, 2655200[':", e.id
	
class QMax_gt_filter(SNPFilter) :
	
	def __init__(self, threshold) :
		self.threshold = threshold
		
	def filter(self, chromosome, dummySRY) :
		if dummySRY.Qmax_gt > self.threshold :
			#other possibilities of return are SequenceInsert(<bases>), SequenceDelete(<length>)
			return SequenceSNP(dummySRY.alt)
		return None #None means keep the reference allele

if __name__ == "__main__" :
	refGenome = Genome(name = 'GRCh37.75_Y-Only')
	persGenome = Genome(name = 'GRCh37.75_Y-Only', SNPs = 'dummySRY', SNPFilter = QMax_gt_filter(10))
	
	geneRef = refGenome.get(Gene, name = 'SRY')[0]
	persGene = persGenome.get(Gene, name = 'SRY')[0]
	
	printing(geneRef)
	print "\n"
	printVs(geneRef, persGene)
	print "\n"
	fancyExonQuery(geneRef)
