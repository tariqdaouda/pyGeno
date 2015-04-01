import time, types, string
import configuration as conf
from rabaDB.rabaSetup import *
from rabaDB.Raba import *
from rabaDB.filters import RabaQuery

def nosave() :
	raise ValueError('You can only save object that are linked to reference genomes')

class pyGenoRabaObject(Raba) :
	"""pyGeno uses rabaDB to persistenly store data. Most persistent 
	objects have classes that inherit from this one (Genome_Raba, 
	Chromosome_Raba, Gene_Raba, Protein_Raba, Exon_Raba). Theses classes 
	are not mean to be accessed directly. Users manipulate wrappers 
	such as : Genome, Chromosome etc... pyGenoRabaObject extends 
	the Raba class by adding a function _curate that is called just 
	before saving. This class is to be considered abstract, and is not 
	meant to be instanciated"""

	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	_raba_abstract = True # not saved in db by default

	def __init__(self) :
		if self is pyGenoRabaObject :
			raise TypeError("This class is abstract")
		
	def _curate(self) :
		"Last operations performed before saving, must be implemented in child"
		raise TypeError("This method is abstract and should be implemented in child")

	def save(self) :
		"""Calls _curate() before performing a normal rabaDB lazy save 
		(saving only occurs if the object has been modified)"""
		
		if self.mutated() :
			self._curate()
		Raba.save(self)

class pyGenoRabaObjectWrapper_metaclass(type) :
	"""This metaclass keeps track of the relationship between wrapped 
	classes and wrappers """
	_wrappers = {}

	def __new__(cls, name, bases, dct) :
		clsObj = type.__new__(cls, name, bases, dct)
		cls._wrappers[dct['_wrapped_class']] = clsObj
		return clsObj

class RLWrapper(object) :
	"""A wrapper for rabalists that replaces raba objects by pyGeno Object"""
	def __init__(self, rabaObj, listObjectType, rl) :
		self.rabaObj = rabaObj
		self.rl = rl
		self.listObjectType = listObjectType

	def __getitem__(self, i) :
		return self.listObjectType(wrapped_object_and_bag = (self.rl[i], self.rabaObj.bagKey))
	
	def __getattr__(self, name) :
		rl =  object.__getattribute__(self, 'rl')
		return getattr(rl, name)

class pyGenoRabaObjectWrapper(object) :
	"""All the wrapper classes such as Genome and Chromosome inherit 
	from this class. It has most that make pyGeno useful, such as 
	get(), count(), ensureIndex(). This class is to be considered 
	abstract, and is not meant to be instanciated"""
	__metaclass__ = pyGenoRabaObjectWrapper_metaclass

	_wrapped_class = None
	_bags = {}

	def __init__(self, wrapped_object_and_bag = (), *args, **kwargs) :
		if self is pyGenoRabaObjectWrapper :
			raise TypeError("This class is abstract")

		if wrapped_object_and_bag != () :
			assert wrapped_object_and_bag[0]._rabaClass is self._wrapped_class
			self.wrapped_object = wrapped_object_and_bag[0]
			self.bagKey = wrapped_object_and_bag[1]
			pyGenoRabaObjectWrapper._bags[self.bagKey][self._getObjBagKey(self.wrapped_object)] = self
		else :
			self.wrapped_object = self._wrapped_class(*args, **kwargs)
			self.bagKey = time.time()
			pyGenoRabaObjectWrapper._bags[self.bagKey] = {}
			pyGenoRabaObjectWrapper._bags[self.bagKey][self._getObjBagKey(self.wrapped_object)] = self

		self._load_sequencesTriggers = set()
		self.loadSequences = True
		self.loadBinarySequences = True

	def _getObjBagKey(self, obj) :
		"""pyGeno objects are kept in bags to ensure that reference 
		objects are loaded only once. This function returns the bag key 
		of the current object"""
		return (obj._rabaClass.__name__, obj.raba_id)

	def _makeLoadQuery(self, objectType, *args, **coolArgs) :
		f = RabaQuery(objectType._wrapped_class, namespace = self._wrapped_class._raba_namespace)
		coolArgs[self._wrapped_class.__name__[:-5]] = self.wrapped_object #[:-5] removes _Raba from class name

		if len(args) > 0 and type(args[0]) is types.ListType :
			for a in args[0] :
				if type(a) is types.DictType :
					f.addFilter(**a)
		else :
			f.addFilter(*args, **coolArgs)

		return f

	def count(self, objectType, *args, **coolArgs) :
		"""Returns the number of elements satisfying the query"""
		return self._makeLoadQuery(objectType, *args, **coolArgs).count()

	def get(self, objectType, *args, **coolArgs) :
		"""Raba Magic inside. This is th function that you use for 
		querying pyGeno's DB.
		
		Usage examples:
		
			* myGenome.get("Gene", name = 'TPST2')
		
			* myGene.get(Protein, id = 'ENSID...')
		
			* myGenome.get(Transcript, {'start >' : x, 'end <' : y})"""
		
		# conf.db.enableDebug(True)
		ret = []
		for e in self._makeLoadQuery(objectType, *args, **coolArgs).iterRun() :
			if issubclass(objectType, pyGenoRabaObjectWrapper) :
				ret.append(objectType(wrapped_object_and_bag = (e, self.bagKey)))
			else :
				ret.append(e)

		return ret

	def iterGet(self, objectType, *args, **coolArgs) :
		"""Same as get. But retuns the elements one by one, much more efficient for large outputs"""

		for e in self._makeLoadQuery(objectType, *args, **coolArgs).iterRun() :
			if issubclass(objectType, pyGenoRabaObjectWrapper) :
				yield objectType(wrapped_object_and_bag = (e, self.bagKey))
			else :
				yield e

	#~ def ensureIndex(self, fields) :
		#~ """
		#~ Warning: May not work on some systems, see ensureGlobalIndex
		#~ 
		#~ Creates a partial index on self (if it does not exist). 
		#~ Ex: myTranscript.ensureIndex('name')"""
		#~ 
		#~ where, whereValues = '%s=?' %(self._wrapped_class.__name__[:-5]), self.wrapped_object
		#~ self._wrapped_class.ensureIndex(fields, where, (whereValues,))

	#~ def dropIndex(self, fields) :
		#~ """Warning: May not work on some systems, see dropGlobalIndex
		#~ 
		#~ Drops a partial index on self. Ex: myTranscript.dropIndex('name')"""

		#~ where, whereValues = '%s=?' %(self._wrapped_class.__name__[:-5]), self.wrapped_object
		#~ self._wrapped_class.dropIndex(fields, where, (whereValues,))
	
	def __getattr__(self, name) :
		"""If a wrapper does not have a specific field, pyGeno will 
		look for it in the wrapped_object"""
		
		if name == 'save' or name == 'delete' :
			raise AttributeError("You can't delete or save an object from wrapper, try .wrapped_object.delete()/save()")
		
		if name in self._load_sequencesTriggers and self.loadSequences :
			self.loadSequences = False
			self._load_sequences()
			return getattr(self, name)
			
		if name[:4] == 'bin_' and self.loadBinarySequences :
			self.updateBinarySequence = False
			self._load_bin_sequence()
			return getattr(self, name)
			
		attr = getattr(self.wrapped_object, name)
		if isRabaObject(attr) :
			attrKey = self._getObjBagKey(attr)
			if attrKey in pyGenoRabaObjectWrapper._bags[self.bagKey] :
				retObj = pyGenoRabaObjectWrapper._bags[self.bagKey][attrKey]
			else :
				wCls = pyGenoRabaObjectWrapper_metaclass._wrappers[attr._rabaClass]
				retObj = wCls(wrapped_object_and_bag = (attr, self.bagKey))
			return retObj
		return attr

	@classmethod
	def getIndexes(cls) :
		"""Returns a list of indexes attached to the object's class. Ex 
		Transcript.getIndexes()"""
		return cls._wrapped_class.getIndexes()

	@classmethod
	def flushIndexes(cls) :
		"""Drops all the indexes attached to the object's class. Ex 
		Transcript.flushIndexes()"""
		return cls._wrapped_class.flushIndexes()
	
	@classmethod
	def help(cls) :
		"""Returns a list of available field for queries. Ex 
		Transcript.help()"""
		return cls._wrapped_class.help().replace("_Raba", "")

	@classmethod
	def ensureGlobalIndex(cls, fields) :
		"""Add a GLOBAL index to the db to speedup lookouts. Fields can be a 
		list of fields for Multi-Column Indices or simply the name of a 
		single field. A global index is an index on the entire type.
		A global index on 'Transcript' on field 'name', will index the names for all the transcripts in the database"""
		cls._wrapped_class.ensureIndex(fields)

	@classmethod
	def dropGlobalIndex(cls, fields) :
		"""Drops an index, the opposite of ensureGlobalIndex()"""
		cls._wrapped_class.dropIndex(fields)

	def _load_sequences(self) :
		"""This lazy abstract function is only called if the object 
		sequences need to be loaded"""
		raise NotImplementedError("This fct loads non binary sequences and should be implemented in child if needed")
	
	def _load_bin_sequence(self) :
		"""Same as _load_sequences(), but loads binary sequences"""
		raise NotImplementedError("This fct loads binary sequences and should be implemented in child if needed")
