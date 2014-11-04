import types

import configuration as conf
import configuration as conf

from pyGenoObjectBases import *
import rabaDB.fields as rf

from tools import UsefulFunctions as uf

"""General guidelines for SNP classes:
----
-All classes must inherit from SNP_INDEL
-All classes name must end with SNP
-A SNP_INDELs must have at least chromosomeNumber, setName, specie, start and ref filled in order to be inserted into sequences
-User can define an alias for the alt field (snp_indel alleles) to indicate the default field from wich to extract alleles
"""

class SequenceSNP_INDEL(object) :
	"This type that a of objects that a SNP filter should return"
	
	class SNPType :
		pass
	class DeletionType :
		pass
	class InsertionType :
		pass

	strTypes = {'SNP' : SNPType, 'DELETION' : DeletionType, 'INSERTION' : InsertionType}
	
	def __init__(self, alleles, polyType, length) :
		"""
		* alleles is string of nucleotides
		* polyType can either be a string among ['SNP', 'DELETION', 'INSERTION'] or one of the types
		SequenceSNP_INDEL.SNPType, SequenceSNP_INDEL.DeletionType, SequenceSNP_INDEL.InsertionType
		* length is the number of nucleatides affected by the ppolymorphism. It must be an int.
		== 1 if it's a SNP, > 1 if its an indel. In that last case it's the number of nucleotides inserted or deleted
		"""
		assert type(length) is types.IntType
		assert length > 0

		if type(alleles) is types.UnicodeType :
			alleles = str(alleles)
		
		if type(polyType) is types.StringType :
			try :
				self.type = self.strTypes[polyType.upper()]
			except KeyError :
				raise TypeError("type, if it's a String, must one of those: ['SNP', 'DELETION', 'INSERTION'], i don't care about the case. You've provided: %s" % polyType)
		elif polyType is SequenceSNP_INDEL.SNPType or polyType is SequenceSNP_INDEL.DeletionType or polyType is SequenceSNP_INDEL.InsertionType :
			self.type = polyType
		else :
			raise TypeError("type, if it's a class, must one of those: SequenceSNP_INDEL.SNPType, SequenceSNP_INDEL.DeletionType, SequenceSNP_INDEL.InsertionType. You've provided: %s" % polyType)
		
		self.alleles = alleles
		self.polyType = polyType
		self.length = length
		self.SNPSources = {}
		
	def addSources(self, sources) :
		"Optional, you can keep a dict that records the polynorphims that were mixed together to make self. They are stored into self.SNPSources"
		self.SNPSources = sources

def defaultSNPFilter(chromosome, **kwargs) :
	"""Default function for filtering snp, does not filter anything. Doesn't apply indels.
	This is also a template that you can use for own filters.
	
	The arguments that will be passed to the function are all polymorphisms that appear at the same position and have the following format:
	SNP_setName1 = snp1, SNP_setName2 = snp2, ...
	
	This allowes you to make custom filters that take into account the snp set.
	
	returns a SequenceSNP_INDEL"""
	
	warn = 'Warning: the default snp filter ignores indels. IGNORED :s of SNP set: %s at pos: %s of chromosome: %s'
	
	alleles = []
	sources = {}
	for snpSet, snp in kwargs.iteritems() :
		pos = snp.start
		if snp.alt[0] == '-' :
			print warn % ('DELETION', snpSet, snp.pos, snp.chromosomeNumber)
		if snp.ref[0] == '-' :
			print warn % ('INSERTION', snpSet, snp.pos, snp.chromosomeNumber)
		else :
			sources[snpSet] = snp
			alleles.append(snp.alt) #if not an indel append the polymorphism
		
	#appends the refence allele to the lot
	refAllele = chromosome.refSequence[pos]
	alleles.append(refAllele)
	
	alleles = uf.encodePolymorphicNucleotide(alleles) #encodes all the polymorphism in a single character
	
	#creates the SequenceSNP_INDEL, length means that it's a snp and not instertion (> 1) or a deletion (< 0)
	ret = SequenceSNP_INDEL(alleles = alleles, polyType = SequenceSNP_INDEL.SNPType, length = 1)
	#optional we keep a record of the polymorphisms that were used during the process
	ret.addSources(sources)
	#print '----', ret.alleles, refAllele, kwargs
	return ret

class SNPMaster(Raba) :
	'This object keeps track of set types'
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	specie = rf.Primitive()
	SNPType = rf.Primitive()
	setName = rf.Primitive()
	_raba_uniques = [('setName',)]

	def _curate(self) :
		self.specie = self.specie.lower()
		self.setName = self.setName.lower()

class SNP_INDEL(pyGenoRabaObject) :
	"All SNPs should inherit from me. The name of the class must end with SNP"
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	_raba_abstract = True # not saved in db

	specie = rf.Primitive()
	setName = rf.Primitive()
	chromosomeNumber = rf.Primitive()

	start = rf.Primitive()
	end = rf.Primitive()
	
	ref = rf.Primitive()
	
	#every SNP_INDEL must have a field alt. This variable allows you to set an alias for it. Chamging the alias does not impact the schema
	altAlias = 'alt'
	
	def __getattribute__(self, k) :
		if k == 'alt' :
			cls = pyGenoRabaObject.__getattribute__(self, '__class__')
			return pyGenoRabaObject.__getattribute__(self, cls.altAlias)
		
		return pyGenoRabaObject.__getattribute__(self, k)

	def __setattr__(self, k, v) :
		if k == 'alt' :
			cls = pyGenoRabaObject.__getattribute__(self, '__class__')
			pyGenoRabaObject.__setattr__(self, cls.altAlias, v)
			
		pyGenoRabaObject.__setattr__(self, k, v)
	
	def _curate(self) :
		self.specie = self.specie.lower()

	@classmethod
	def ensureGlobalIndex(cls, fields) :
		cls.ensureIndex(fields)
	
class CasavaSNP(SNP_INDEL) :
	"A SNP of Casava"
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	
	alleles = rf.Primitive()
	bcalls_used = rf.Primitive()
	bcalls_filt = rf.Primitive()
	QSNP = rf.Primitive()
	Qmax_gt = rf.Primitive()
	max_gt_poly_site = rf.Primitive()
	Qmax_gt_poly_site = rf.Primitive()
	A_used = rf.Primitive()
	C_used = rf.Primitive()
	G_used = rf.Primitive()
	T_used = rf.Primitive()
	
	altAlias = 'alleles'

class dbSNPSNP(SNP_INDEL) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	
	# To add/remove a field comment/uncomentd it. Beware, adding or removing a field results in a significant overhead the first time you relaunch pyGeno. You may have to delete and reimport some snps sets.
	
	#####RSPOS = rf.Primitive() #Chr position reported in already saved into field start
	RS = rf.Primitive() #dbSNP ID (i.e. rs number)
	ALT =  rf.Primitive()
	RV = rf.Primitive() #RS orientation is reversed
	#VP = rf.Primitive() #Variation Property.  Documentation is at ftp://ftp.ncbi.nlm.nih.gov/snp/specs/dbSNP_BitField_latest.pdf
	#GENEINFO = rf.Primitive() #Pairs each of gene symbol:gene id.  The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)
	dbSNPBuildID = rf.Primitive() #First dbSNP Build for RS
	#SAO = rf.Primitive() #Variant Allele Origin: 0 - unspecified, 1 - Germline, 2 - Somatic, 3 - Both
	#SSR = rf.Primitive() #Variant Suspect Reason Codes (may be more than one value added together) 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other
	#WGT = rf.Primitive() #Weight, 00 - unmapped, 1 - weight 1, 2 - weight 2, 3 - weight 3 or more
	VC = rf.Primitive() #Variation Class
	PM = rf.Primitive() #Variant is Precious(Clinical,Pubmed Cited)
	#TPA = rf.Primitive() #Provisional Third Party Annotation(TPA) (currently rs from PHARMGKB who will give phenotype data)
	#PMC = rf.Primitive() #Links exist to PubMed Central article
	#S3D = rf.Primitive() #Has 3D structure - SNP3D table
	#SLO = rf.Primitive() #Has SubmitterLinkOut - From SNP->SubSNP->Batch.link_out
	#NSF = rf.Primitive() #Has non-synonymous frameshift A coding region variation where one allele in the set changes all downstream amino acids. FxnClass = 44
	#NSM = rf.Primitive() #Has non-synonymous missense A coding region variation where one allele in the set changes protein peptide. FxnClass = 42
	#NSN = rf.Primitive() #Has non-synonymous nonsense A coding region variation where one allele in the set changes to STOP codon (TER). FxnClass = 41
	#REF = rf.Primitive() #Has reference A coding region variation where one allele in the set is identical to the reference sequence. FxnCode = 8
	#SYN = rf.Primitive() #Has synonymous A coding region variation where one allele in the set does not change the encoded amino acid. FxnCode = 3
	#U3 = rf.Primitive() #In 3' UTR Location is in an untranslated region (UTR). FxnCode = 53
	#U5 = rf.Primitive() #In 5' UTR Location is in an untranslated region (UTR). FxnCode = 55
	#ASS = rf.Primitive() #In acceptor splice site FxnCode = 73
	#DSS = rf.Primitive() #In donor splice-site FxnCode = 75
	#INT = rf.Primitive() #In Intron FxnCode = 6
	#R3 = rf.Primitive() #In 3' gene region FxnCode = 13
	#R5 = rf.Primitive() #In 5' gene region FxnCode = 15
	#OTH = rf.Primitive() #Has other variant with exactly the same set of mapped positions on NCBI refernce assembly.
	#CFL = rf.Primitive() #Has Assembly conflict. This is for weight 1 and 2 variant that maps to different chromosomes on different assemblies.
	#ASP = rf.Primitive() #Is Assembly specific. This is set if the variant only maps to one assembly
	MUT = rf.Primitive() #Is mutation (journal citation, explicit fact): a low frequency variation that is cited in journal and other reputable sources
	VLD = rf.Primitive() #Is Validated.  This bit is set if the variant has 2+ minor allele count based on frequency or genotype data.
	G5A = rf.Primitive() #>5% minor allele frequency in each and all populations
	G5 = rf.Primitive() #>5% minor allele frequency in 1+ populations
	#HD = rf.Primitive() #Marker is on high density genotyping kit (50K density or greater).  The variant may have phenotype associations present in dbGaP.
	#GNO = rf.Primitive() #Genotypes available. The variant has individual genotype (in SubInd table).
	KGValidated = rf.Primitive() #1000 Genome validated
	#KGPhase1 = rf.Primitive() #1000 Genome phase 1 (incl. June Interim phase 1)
	#KGPilot123 = rf.Primitive() #1000 Genome discovery all pilots 2010(1,2,3)
	#KGPROD = rf.Primitive() #Has 1000 Genome submission
	OTHERKG = rf.Primitive() #non-1000 Genome submission
	#PH3 = rf.Primitive() #HAP_MAP Phase 3 genotyped: filtered, non-redundant
	#CDA = rf.Primitive() #Variation is interrogated in a clinical diagnostic assay
	#LSD = rf.Primitive() #Submitted from a locus-specific database
	#MTP = rf.Primitive() #Microattribution/third-party annotation(TPA:GWAS,PAGE)
	#OM = rf.Primitive() #Has OMIM/OMIA
	#NOC = rf.Primitive() #Contig allele not present in variant allele list. The reference sequence allele at the mapped position is not present in the variant allele list, adjusted for orientation.
	#WTD = rf.Primitive() #Is Withdrawn by submitter If one member ss is withdrawn by submitter, then this bit is set.  If all member ss' are withdrawn, then the rs is deleted to SNPHistory
	#NOV = rf.Primitive() #Rs cluster has non-overlapping allele sets. True when rs set has more than 2 alleles from different submissions and these sets share no alleles in common.
	#CAF = rf.Primitive() #An ordered, comma delimited list of allele frequencies based on 1000Genomes, starting with the reference allele followed by alternate alleles as ordered in the ALT column. Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previuosly reported in VCF as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter
	COMMON = rf.Primitive() #RS is a common SNP.  A common SNP is one that has at least one 1000Genomes population with a minor allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.
	#CLNHGVS = rf.Primitive() #Variant names from HGVS.    The order of these variants corresponds to the order of the info in the other clinical  INFO tags.
	#CLNALLE = rf.Primitive() #Variant alleles from REF or ALT columns.  0 is REF, 1 is the first ALT allele, etc.  This is used to match alleles with other corresponding clinical (CLN) INFO tags.  A value of -1 indicates that no allele was found to match a corresponding HGVS allele name.
	#CLNSRC = rf.Primitive() #Variant Clinical Chanels
	#CLNORIGIN = rf.Primitive() #Allele Origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other
	#CLNSRCID = rf.Primitive() #Variant Clinical Channel IDs
	#CLNSIG = rf.Primitive() #Variant Clinical Significance, 0 - unknown, 1 - untested, 2 - non-pathogenic, 3 - probable-non-pathogenic, 4 - probable-pathogenic, 5 - pathogenic, 6 - drug-response, 7 - histocompatibility, 255 - other
	#CLNDSDB = rf.Primitive() #Variant disease database name
	#CLNDSDBID = rf.Primitive() #Variant disease database ID
	#CLNDBN = rf.Primitive() #Variant disease name
	#CLNACC = rf.Primitive() #Variant Accession and Versions
	
	#FILTER_NC = rf.Primitive() #Inconsistent Genotype Submission For At Least One Sample
	
	altAlias = 'ALT'

class TopHatSNP(SNP_INDEL) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	pass
