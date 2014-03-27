import time, types, string
import configuration as conf

from rabaDB.setup import *
from rabaDB.Raba import *
from rabaDB.filters import RabaQuery

def nosave() :
	raise ValueError('You can only save object that are linked to reference genomes')

class pyGenoRabaObject(Raba) :
	_raba_namespace = conf.pyGeno_RABA_NAMESPACE
	_raba_abstract = True # not saved in db

	def __init__(self) :
		if self is pyGenoRabaObject :
			raise TypeError("This class is abstract")

	def _curate(self) :
		"last operations performed before saving, must be implemented in child"
		raise TypeError("This method is abstract and should be implemented in child")

	def save(self) :
		if self.mutated() :
			self._curate()
		Raba.save(self)

class pyGenoRabaObjectWrapper_metaclass(type) :

	_wrappers = {}

	def __new__(cls, name, bases, dct) :
		clsObj = type.__new__(cls, name, bases, dct)
		cls._wrappers[dct['_wrapped_class']] = clsObj
		return clsObj

class pyGenoRabaObjectWrapper(object) :

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

		self.loadBinarySequences = True

	def _getObjBagKey(self, obj) :
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

	def get(self, objectType, *args, **coolArgs) :
		"""Raba Magic inside. Loads anything with any parameters. Don't worry about memory Raba takes care of it
		usage:
		load("Gene", name = 'TPST2')
		load(Gene, name = 'TPST2')
		load(Transcript, {'len >' : 60})

		You can use string as the object type to avoid circular imports"""
		ret = []
		for e in self._makeLoadQuery(objectType, *args, **coolArgs).iterRun() :
			ret.append(objectType(wrapped_object_and_bag = (e, self.bagKey)))
		return ret

	def iterGet(self, objectType, *args, **coolArgs) :
		"""Same as load. But retuns the elements one by one, much more efficient for large outputs"""

		for e in self._makeLoadQuery(objectType, *args, **coolArgs).iterRun() :
			yield objectType(wrapped_object_and_bag = (e, self.bagKey))

	def ensureIndex(self, objectType, *fields) :
		where, whereValues = '%s=?' %(self._wrapped_class.__name__[:-5]), self.wrapped_object
		objectType._wrapped_class.ensureIndex(fields, where, (whereValues,))

	def __getattr__(self, name) :
		if name == 'save' or name == 'delete' :
			raise AttributeError("You can't delete or save an object from wrapper, try .wrapped_object.delete()/save()")

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
	def help(cls) :
		return cls._wrapped_class.help()

	@classmethod
	def enureGlobalIndex(cls, fields) :
		"""Add a GLOBAL index to the db to speedup lookouts. Fields can be a list of fields for Multi-Column Indices or simply the name of a single field
		A global index is an index on the entire type. Ex global index on transcript field name, will index the names for all the transcripts in the database"""
		cls._wrapped_class.enureIndex(fields)

	@classmethod
	def dropGlobalIndex(cls, fields) :
		cls._wrapped_class.dropIndex(fields)

	def _load_bin_sequence(self) :
		raise NotImplementedError("This fct loads binary sequences and should be implemented in child")

	def __getattribute__(self, name) :
		if name[:4] == 'bin_' and self.loadBinarySequences :
			self._load_bin_sequence()
			self.updateBinarySequence = False

		return object.__getattribute__(self, name)

