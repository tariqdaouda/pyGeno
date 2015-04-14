import types

import configuration as conf

from pyGenoObjectBases import *
import rabaDB.fields as rf

# from tools import UsefulFunctions as uf

"""General guidelines for SNP classes:
----
-All classes must inherit from SNP_INDEL
-All classes name must end with SNP
-A SNP_INDELs must have at least chromosomeNumber, setName, species, start and ref filled in order to be inserted into sequences
-User can define an alias for the alt field (snp_indel alleles) to indicate the default field from wich to extract alleles
"""

def getSNPSetsList() :
	"""Return the names of all imported snp sets"""
	import rabaDB.filters as rfilt
	f = rfilt.RabaQuery(SNPMaster)
	names = []
	for g in f.iterRun() :
		names.append(g.setName)
	return names

class SNPMaster(Raba) :
	'This object keeps track of SNP sets and their types'
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	species = rf.Primitive()
	SNPType = rf.Primitive()
	setName = rf.Primitive()
	_raba_uniques = [('setName',)]

	def _curate(self) :
		self.species = self.species.lower()
		self.setName = self.setName.lower()

class SNP_INDEL(Raba) :
	"All SNPs should inherit from me. The name of the class must end with SNP"
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	_raba_abstract = True # not saved in db

	species = rf.Primitive()
	setName = rf.Primitive()
	chromosomeNumber = rf.Primitive()

	start = rf.Primitive()
	end = rf.Primitive()
	
	ref = rf.Primitive()
	
	#every SNP_INDEL must have a field alt. This variable allows you to set an alias for it. Chamging the alias does not impact the schema
	altAlias = 'alt'
	
	def __init__(self, *args, **kwargs) :
		pass

	def __getattribute__(self, k) :
		if k == 'alt' :
			cls = Raba.__getattribute__(self, '__class__')
			return Raba.__getattribute__(self, cls.altAlias)
		
		return Raba.__getattribute__(self, k)

	def __setattr__(self, k, v) :
		if k == 'alt' :
			cls = Raba.__getattribute__(self, '__class__')
			Raba.__setattr__(self, cls.altAlias, v)
			
		Raba.__setattr__(self, k, v)
	
	def _curate(self) :
		self.species = self.species.lower()

	@classmethod
	def ensureGlobalIndex(cls, fields) :
		cls.ensureIndex(fields)

	def __repr__(self) :
		return "%s> chr: %s, start: %s, end: %s, alt: %s, ref: %s" %(self.__class__.__name__, self.chromosomeNumber, self.start, self.end, self.alleles, self.ref)
	
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

class AgnosticSNP(SNP_INDEL) :
	"""This is a generic SNPs/Indels format that you can easily make from the result of any SNP caller. AgnosticSNP files are tab delimited files such as:

	chromosomeNumber	uniqueId  start	      end	   ref    alleles	quality	 caller
	Y					   1 	 2655643	2655644		T		AG	     30		 TopHat
	Y					   2 	 2655645	2655647		-		AG	     28		 TopHat
	Y					   3 	 2655648	2655650		TT		-	     10		 TopHat

	All positions must be 0 based
	The '-' indicates a deletion or an insertion. Collumn order has no importance.
	"""

	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	
	alleles = rf.Primitive()
	quality = rf.Primitive()
	caller = rf.Primitive()
	uniqueId = rf.Primitive() # polymorphism id
	
	altAlias = 'alleles'

	def __repr__(self) :
		return "AgnosticSNP> start: %s, end: %s, quality: %s, caller %s, alt: %s, ref: %s" %(self.start, self.end, self.quality, self.caller, self.alleles, self.ref)

class dbSNPSNP(SNP_INDEL) :
	"This class is for SNPs from dbSNP. Feel free to uncomment the fields that you need"
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
	"A SNP from Top Hat, not implemented"
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	pass
