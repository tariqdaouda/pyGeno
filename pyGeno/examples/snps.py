from pyGeno.Genome import Genome
from pyGeno.Gene import Gene
from pyGeno.Transcript import Transcript
from pyGeno.Protein import Protein
from pyGeno.SNP import SequenceSNP_INDEL

def mySNPFilter(chromosome, dummySRY) :	
	warn = 'Warning: the default snp filter ignores indels. IGNORED %s of SNP set: %s at pos: %s of chromosome: %s'
	#~ print "aaaaaaaaaaaaaaaaaaaaaaa", dummySRY.alt, "pppp"
	
	if dummySRY.alt == "R" :
		#~ ret = SequenceSNP_INDEL(alleles = alleles, polyType = SequenceSNP_INDEL.SNPType, length = 1)
		#~ return SequenceSNP_INDEL(alleles = '', polyType = SequenceSNP_INDEL.DeletionType, length = 3)
		return SequenceSNP_INDEL(alleles = 'ATG', polyType = SequenceSNP_INDEL.InsertionType, length = 1)
		
	#~ refAllele = chromosome.refSequence[pos]
	#~ alleles.append(refAllele)
	#~ 
	#~ alleles = uf.encodePolymorphicNucleotide(alleles) #encodes all the polymorphism in a single character
	#~ 
	#~ 
	#~ ret.addSources(sources)
	#~ return ret

ref = Genome(name = 'GRCh37.75_Y-Only')
dummy = Genome(name = 'GRCh37.75_Y-Only', SNPs = 'dummySRY', SNPFilter = mySNPFilter)
geneRef = ref.get(Gene, name = 'SRY')[0]

#~ print "printing reference sequences\n-------"
#~ for trans in geneRef.get(Transcript) :
	#~ print "\t-Transcript name:", trans.name
	#~ print "\t-Protein:", trans.protein.sequence
	#~ print "\t-Exons:"
	#~ for e in trans.exons :
		#~ print "\t\t Number:", e.number
		#~ print "\t\t-CDS:", e.CDS
		#~ print "\t\t-Strand:", e.strand
		#~ print "\t\t-CDS_start:", e.CDS_start
		#~ print "\t\t-CDS_end:", e.CDS_end
#~ 
#~ print '\n======\n'

print "Vs personalized sequences\n------"

geneDummy = dummy.get(Gene, name = 'SRY')[0]
for trans in geneDummy.get(Transcript) :
	refProt = ref.get(Protein, id = trans.protein.id)[0]
	print trans.protein.id
	print "\tref:" + refProt.sequence[:20] + "..."
	print "\tper:" + trans.protein.sequence[:20] + "..."
	print
