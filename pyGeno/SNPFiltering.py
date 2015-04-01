import types

import configuration as conf

from tools import UsefulFunctions as uf

class Sequence_modifiers(object) :
	"""Abtract Class. All sequence must inherit from me"""
	def __init__(self, sources = {}) :
		self.sources = sources

	def addSource(self, name, snp) :
		"Optional, you can keep a dict that records the polymorphims that were mixed together to make self. They are stored into self.sources"
		self.sources[name] = snp

class SequenceSNP(Sequence_modifiers) :
	"""Represents a SNP to be applied to the sequence"""
	def __init__(self, alleles, sources = {}) :
		Sequence_modifiers.__init__(self, sources)
		if type(alleles) is types.ListType :
			self.alleles = uf.encodePolymorphicNucleotide(''.join(alleles))
		else :
			self.alleles = uf.encodePolymorphicNucleotide(alleles)
	
class SequenceInsert(Sequence_modifiers) :
	"""Represents an Insertion to be applied to the sequence"""
	
	def __init__(self, bases, sources = {}) :
		Sequence_modifiers.__init__(self, sources)
		self.bases = bases

class SequenceDel(Sequence_modifiers) :
	"""Represents a Deletion to be applied to the sequence"""
	def __init__(self, length, sources = {}) :
		Sequence_modifiers.__init__(self, sources)
		self.length = length

class SNPFilter(object) :
	"""Abtract Class. All filters must inherit from me"""
	
	def __init__(self) :
		pass

	def filter(self, chromosome, **kwargs) :
		raise NotImplemented("Must be implemented in child")

class DefaultSNPFilter(SNPFilter) :
	"""
	Default filtering object, does not filter anything. Doesn't apply indels.
	This is also a template that you can use for own filters. A prototype for a custom filter might be::
	
		class MyFilter(SNPFilter) :
			def __init__(self, thres) :
				self.thres = thres
			
			def filter(chromosome, SNP_Set1 = None, SNP_Set2 = None ) :
				if SNP_Set1.alt is not None and (SNP_Set1.alt == SNP_Set2.alt) and SNP_Set1.Qmax_gt > self.thres :
					return SequenceSNP(SNP_Set1.alt)
				return None
			
	Where SNP_Set1 and SNP_Set2 are the actual names of the snp sets supplied to the genome. In the context of the function
	they represent single polymorphisms derived from thoses sets that occur at the same position.

	Whatever goes on into the function is absolutely up to you, but in the end, it must return an object of one of the following classes:

		* SequenceSNP

		* SequenceInsert

		* SequenceDel

		* None

		"""

	def __init__(self) :
		SNPFilter.__init__(self)

	def filter(self, chromosome, **kwargs) :
		"""The default filter mixes applied all SNPs and ignores Insertions and Deletions."""
		warn = 'Warning: the default snp filter ignores indels. IGNORED %s of SNP set: %s at pos: %s of chromosome: %s'
		
		sources = {}
		alleles = []
		for snpSet, snp in kwargs.iteritems() :
			pos = snp.start
			if snp.alt[0] == '-' :
				pass
				# print warn % ('DELETION', snpSet, snp.start, snp.chromosomeNumber)
			elif snp.ref[0] == '-' :
				pass
				# print warn % ('INSERTION', snpSet, snp.start, snp.chromosomeNumber)
			else :
				sources[snpSet] = snp
				alleles.append(snp.alt) #if not an indel append the polymorphism
			
		#appends the refence allele to the lot
		refAllele = chromosome.refSequence[pos]
		alleles.append(refAllele)
		sources['ref'] = refAllele

		#optional we keep a record of the polymorphisms that were used during the process
		return SequenceSNP(alleles, sources = sources)
