import random, copy, types

def aux_insertTree(childTree, parentTree):
	"""This a private (You shouldn't have to call it) recursive function that inserts a child tree into a parent tree."""
	if childTree.x1 != None and childTree.x2 != None :
		parentTree.insert(childTree.x1, childTree.x2, childTree.name, childTree.referedObject)

	for c in childTree.children:
		aux_insertTree(c, parentTree)
		
def aux_moveTree(offset, tree):
	"""This a private recursive (You shouldn't have to call it) function that translates tree(and it's children) to a given x1"""
	if tree.x1 != None and tree.x2 != None :
		tree.x1, tree.x2 = tree.x1+offset, tree.x2+offset
		
	for c in tree.children:
		aux_moveTree(offset, c)
		
class SegmentTree :
	""" Optimised genome annotations.
	A segment tree is an arborescence of segments. First position is inclusive, second exlusive, respectively refered to as x1 and x2.
	A segment tree has the following properties :
	
	* The root has no x1 or x2 (both set to None).
	
	* Segment are arrangend in an ascending order
	
	* For two segment S1 and S2 : [S2.x1, S2.x2[ C [S1.x1, S1.x2[ <=> S2 is a child of S1
	
	Here's an example of a tree :
	
	* Root : 0-15
	
	* ---->Segment : 0-12
	
	* ------->Segment : 1-6
	
	* ---------->Segment : 2-3
	
	* ---------->Segment : 4-5
	
	* ------->Segment : 7-8
	
	* ------->Segment : 9-10
	
	* ---->Segment : 11-14
	
	* ------->Segment : 12-14
	
	* ---->Segment : 13-15
	
	Each segment can have a 'name' and a 'referedObject'. ReferedObject are objects are stored within the graph for future usage.
	These objects are always stored in lists. If referedObject is already a list it will be stored as is.
	"""
	
	def __init__(self, x1 = None, x2 = None, name = '', referedObject = [], father = None, level = 0) :
		if x1 > x2 :
			self.x1, self.x2 = x2, x1
		else :
			self.x1, self.x2 = x1, x2
		
		self.father = father
		self.level = level
		self.id = random.randint(0, 10**8)
		self.name = name
		
		self.children = []
		self.referedObject = referedObject

	def __addChild(self, segmentTree, index = -1) :
		segmentTree.level = self.level + 1
		if index < 0 :
			self.children.append(segmentTree)
		else :
			self.children.insert(index, segmentTree)
	
	def insert(self, x1, x2, name = '', referedObject = []) :
		"""Insert the segment in it's right place and returns it. 
		If there's already a segment S as S.x1 == x1 and S.x2 == x2. S.name will be changed to 'S.name U name' and the
		referedObject will be appended to the already existing list"""
		
		if x1 > x2 :
			xx1, xx2 = x2, x1
		else :
			xx1, xx2 = x1, x2

		rt = None
		insertId = None
		childrenToRemove = []
		for i in range(len(self.children)) :
			if self.children[i].x1 == xx1 and xx2 == self.children[i].x2 :
				self.children[i].name = self.children[i].name + ' U ' + name
				self.children[i].referedObject.append(referedObject)
				return self.children[i]
			
			if self.children[i].x1 <= xx1 and xx2 <= self.children[i].x2 :
				return self.children[i].insert(x1, x2, name, referedObject)
			
			elif xx1 <= self.children[i].x1 and self.children[i].x2 <= xx2 :
				if rt == None :
					if type(referedObject) is types.ListType :
						rt = SegmentTree(xx1, xx2, name, referedObject, self, self.level+1)
					else :
						rt = SegmentTree(xx1, xx2, name, [referedObject], self, self.level+1)
					
					insertId = i
					
				rt.__addChild(self.children[i])
				self.children[i].father = rt
				childrenToRemove.append(self.children[i])
		
			elif xx1 <= self.children[i].x1 and xx2 <= self.children[i].x2 :
				insertId = i
				break
				
		if rt != None :
			self.__addChild(rt, insertId)
			for c in childrenToRemove :
				self.children.remove(c)
		else :
			if type(referedObject) is types.ListType :
				rt = SegmentTree(xx1, xx2, name, referedObject, self, self.level+1)
			else :
				rt = SegmentTree(xx1, xx2, name, [referedObject], self, self.level+1)
			
			if insertId != None :
				self.__addChild(rt, insertId)
			else :
				self.__addChild(rt)
		
		return rt
	
	def insertTree(self, childTree):
		"""inserts childTree in the right position (regions will be rearanged to fit the organisation of self)"""
		aux_insertTree(childTree, self)
		
	#~ def included_todo(self, x1, x2=None) :
		#~ "Returns all the segments where [x1, x2] is included"""
		#~ pass
		
	def intersect(self, x1, x2 = None) :
		"""Returns a list of all segments intersected by [x1, x2]"""
		
		def condition(x1, x2, tree) :
			#print self.id, tree.x1, tree.x2, x1, x2
			if (tree.x1 != None and tree.x2 != None) and (tree.x1 <= x1 and x1 < tree.x2 or tree.x1 <= x2 and x2 < tree.x2) :
				return True
			return False
			
		if x2 == None :
			xx1, xx2 = x1, x1
		elif x1 > x2 :
			xx1, xx2 = x2, x1
		else :
			xx1, xx2 = x1, x2
			
		c1 = self.__dichotomicSearch(xx1)
		c2 = self.__dichotomicSearch(xx2)
		
		if c1 == -1 or c2 == -1 :
			return []
			
		if xx1 < self.children[c1].x1 :
			c1 -= 1
			
		inter = self.__radiateDown(x1, x2, c1, condition)
		if self.children[c1].id == self.children[c2].id :
			inter.extend(self.__radiateUp(x1, x2, c2+1, condition))
		else :
			inter.extend(self.__radiateUp(x1, x2, c2, condition))
		
		ret = []
		for c in inter :
			ret.extend(c.intersect(x1, x2))
		
		inter.extend(ret)
		return inter
		
	def __dichotomicSearch(self, x1) :
		r1 = 0
		r2 = len(self.children)-1
		pos = -1
		while (r1 <= r2) :
			pos = (r1+r2)/2
			val = self.children[pos].x1

			if val == x1 :
				return pos
			elif x1 < val :
				r2 = pos -1
			else :
				r1 = pos +1
		
		return pos
	
	def __radiateDown(self, x1, x2, childId, condition) :
		"Radiates down: walks self.children downward until condition is no longer verifed or there's no childrens left "
		ret = []
		i = childId
		while 0 <= i :
			if condition(x1, x2, self.children[i]) :
				ret.append(self.children[i])
			else :
				break
			i -= 1
		return ret
	
	def __radiateUp(self, x1, x2, childId, condition) :
		"Radiates uo: walks self.children upward until condition is no longer verifed or there's no childrens left "
		ret = []
		i = childId
		while i < len(self.children):
			if condition(x1, x2, self.children[i]) :
				ret.append(self.children[i])
			else :
				break
			i += 1
		return ret
	
	def emptyChildren(self) :
		"""Kills of children"""
		self.children = []
	
	def removeGaps(self) :
		"""Remove all gaps between regions"""
		
		for i in range(1, len(self.children)) :
			if self.children[i].x1 > self.children[i-1].x2:				
				aux_moveTree(self.children[i-1].x2-self.children[i].x1, self.children[i])
		
	def getX1(self) :
		"""Returns the starting position of the tree"""
		if self.x1 != None :
			return self.x1
		return self.children[0].x1

	def getX2(self) :
		"""Returns the ending position of the tree"""
		if self.x2 != None :
			return self.x2
		return self.children[-1].x2
	
	def getIndexedLength(self) :
		"""Returns the total length of indexed regions"""
		if self.x1 != None and self.x2 != None:
			return self.x2 - self.x1
		else :
			if len(self.children) == 0 :
				return 0
			else :
				l = self.children[0].x2 - self.children[0].x1
				for i in range(1, len(self.children)) :
					l += self.children[i].x2 - self.children[i].x1 - max(0, self.children[i-1].x2 - self.children[i].x1)
				return l
	
	def getFirstLevel(self) :
		"""returns a list of couples (x1, x2) of all the first level indexed regions"""
		res = []
		if len(self.children) > 0 :
			for c in self.children:
				res.append((c.x1, c.x2)) 
		else :
			if self.x1 != None :
				res = [(self.x1, self.x2)]
			else :
				res = None
		return res
		
	def flatten(self) :
		"""Flattens the tree. The tree become a tree of depth 1 where overlapping regions have been merged together"""
		if len(self.children) > 1 :
			children = self.children
			self.emptyChildren()
			
			children[0].emptyChildren()
			x1 = children[0].x1
			x2 = children[0].x2
			refObjs = [children[0].referedObject]
			name = children[0].name
			
			for i in range(1, len(children)) :
				children[i].emptyChildren()
				if children[i-1] >= children[i] :
					x2 = children[i].x2
					refObjs.append(children[i].referedObject)
					name += " U " + children[i].name
				else :
					if len(refObjs) == 1 :
						refObjs = refObjs[0]
		
					self.insert(x1, x2, name, refObjs)
					x1 = children[i].x1
					x2 = children[i].x2
					refObjs = [children[i].referedObject]
					name = children[i].name
			
			if len(refObjs) == 1 :
				refObjs = refObjs[0]
		
			self.insert(x1, x2, name, refObjs)

	def move(self, newX1) :
		"""Moves tree to a new starting position, updates x1s of children"""
		if self.x1 != None and self.x2 != None :
			offset = newX1-self.x1
			aux_moveTree(offset, self)
		elif len(self.children) > 0 :
			offset = newX1-self.children[0].x1
			aux_moveTree(offset, self)

	def __str__(self) :
		strRes = self.__str()
		
		offset = ''
		for i in range(self.level+1) :
			offset += '\t'
			
		for c in self.children :
			strRes += '\n'+offset+'-->'+str(c)
		
		return strRes
	
	def __str(self) :
		if self.x1 == None :
			if len(self.children) > 0 :
				return "Root: %d-%d, name: %s, id: %d, obj: %s" %(self.children[0].x1, self.children[-1].x2, self.name, self.id, repr(self.referedObject))
			else :
				return "Root: EMPTY , name: %s, id: %d, obj: %s" %(self.name, self.id, repr(self.referedObject))
		else :
			return "Segment: %d-%d, name: %s, id: %d, father id: %d, obj: %s" %(self.x1, self.x2, self.name, self.id, self.father.id, repr(self.referedObject))
			
	
	def __len__(self) :
		"returns the size of the complete indexed region"
		if self.x1 != None and self.x2 != None :
			return self.x2-self.x1
		else :
			return self.children[-1].x2 - self.children[0].x1
	
	def __repr__(self):
		return 'Segment Tree, id:%s, father id:%s, (x1, x2): (%s, %s)' %(self.id, self.father.id, self.x1, self.x2)


if __name__== "__main__" :
	s = SegmentTree()
	s.insert(5, 10, 'region 1')
	s.insert(8, 12, 'region 2')
	s.insert(5, 8, 'region 3')
	s.insert(34, 40, 'region 4')
	s.insert(35, 38, 'region 5')
	s.insert(36, 37, 'region 6', 'aaa')
	s.insert(36, 37, 'region 6', 'aaa2')
	print "Tree:"
	print s
	print "indexed length", s.getIndexedLength()
	print "removing gaps and adding region 7 : [13-37["
	s.removeGaps()
	#s.insert(13, 37, 'region 7')
	print s
	print "indexed length", s.getIndexedLength()
	#print "intersections"
	#for c in [6, 10, 14, 1000] :
	#	print c, s.intersect(c)
	
	print "Move"
	s.move(0)
	print s
	print "indexed length", s.getIndexedLength()
