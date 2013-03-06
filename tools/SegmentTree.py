import random, copy

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
	"""This is like an R-tree but with segments
	Root : 0-15
	---->Segment : 0-12
	------->Segment : 1-6
	---------->Segment : 2-3
	---------->Segment : 4-5
	------->Segment : 7-8
	------->Segment : 9-10
	---->Segment : 11-14
	------->Segment : 12-14
	---->Segment : 13-15
	"""
	
	def __init__(self, x1 = None, x2 = None, name = '', referedObject = None, father = None, level = 0) :
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
	
	def insert(self, x1, x2, name = '', referedObject = None) :
		"""Insert the segment at its correct place and returns it"""
		if x1 > x2 :
			xx1, xx2 = x2, x1
		else :
			xx1, xx2 = x1, x2

		rt = None
		insertId = None
		childrenToRemove = []
		for i in range(len(self.children)) :
			if self.children[i].x1 == xx1 and xx2 == self.children[i].x2 :
				return self.children[i]
			
			if self.children[i].x1 <= xx1 and xx2 <= self.children[i].x2 :
				return self.children[i].insert(x1, x2, name, referedObject)
			
			elif xx1 <= self.children[i].x1 and self.children[i].x2 <= xx2 :
				if rt == None :
					rt = SegmentTree(xx1, xx2, name, referedObject, self, self.level+1)
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
			rt = SegmentTree(xx1, xx2, name, referedObject, self, self.level+1)
			if insertId != None :
				self.__addChild(rt, insertId)
			else :
				self.__addChild(rt)
		
		return rt
	
	def insertTree(self, childTree):
		"""inserts segTree in the right position (regions will rearanged to fit the organisation of self)"""
		aux_insertTree(childTree, self)
	
	def intersect(self, x1, x2= None) :
		"""Returns all the segments intersected by x1, x2"""
		if x2 == None :
			xx1, xx2 = x1, x1
		elif x1 > x2 :
			xx1, xx2 = x2, x1
		else :
			xx1, xx2 = x1, x2
		
		ret = []
		if (self.x1 != None and self.x2 != None) and (self.x1 <= xx1 and xx2 <= self.x2) :
			ret.append(self)

		for c in self.children :
			ret.extend(c.intersect(x1, x2))

		return ret

	def emptyChildren(self) :
		self.children = []
		
	def removeGaps_bck(self) :
		"""returns a tree similar to self but where the gaps between succesive indexed regions have been removed:
		Regions are moved so they adajcent to the one before, the transfornations also affect children"""
		res = SegmentTree()
		
		cc = copy.copy(self.children[0])
		res.insertTree(cc)
		previousC = cc
		for i in range(1, len(self.children)) :
			cc = copy.copy(self.children[i])
			if self.children[i].x1 > previousC.x2:				
				aux_moveTree(previousC.x2-cc.x1, cc)
			res.insertTree(cc)
			previousC = cc
		return res
	
	def removeGaps(self) :
		"""Remove all gaps between regions"""
		
		for i in range(1, len(self.children)) :
			if self.children[i].x1 > self.children[i-1].x2:				
				aux_moveTree(self.children[i-1].x2-self.children[i].x1, self.children[i])
		
	def getX1(self) :
		if self.x1 != None :
			return self.x1
		return self.children[0].x1

	def getX2(self) :
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
		'returns a list of couples (x1, x2) of all the first level indexed regions'
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


	def flattened_bck(self) :
		r"""Returns a flattend version of the tree. A list of SegmentTrees where all overlapping children have been merged together"""
		root = SegmentTree()
		if self.x1 == None and self.x2 == None:
			if len(self.children) == 1 :
				root.insert(self.children[0].x1, self.children[0].x2, '', [self.children[0].referedObject])
			elif len(self.children) > 1 :
				x1 = self.children[0].x1
				x2 = self.children[0].x2
				refObjs = [self.children[0].referedObject]
				for c in self.children[1:]:
					if x2 >= c.x1 :
						x2 = c.x2
						refObjs.append(c.referedObject)
					else :
						root.insert(x1, x2, '', referedObject = refObjs)
						x1 = c.x1
						x2 = c.x2
						refObjs = [c.referedObject]
						
				root.insert(x1, x2, '', referedObject = refObjs)
		else :
			return self
		return root
		
	def flatten(self) :
		r"""Flattens the tree. The tree become a tree of depth 1 where overlapping regions have been merged together"""
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
			aux_moveTree(newX1-self.x1, self)
		elif len(self.children) > 0:
			offset = newX1-self.children[0].x1
			for c in self.children :
				aux_moveTree(offset, c)

	def translate(self, offset) :
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
		return self.getIndexedLength()
		
	def __len__bck(self) :
		r"""if root returns whole region of interest length. if not it is identical to getIndexedRegionsLength()"""
		
		if self.x1 != None :
			return self.getIndexedRegionsLength()
		else :
			if len(self.children) == 0 :
				return 0
			else :
				xx1, xx2 = self.children[0].x1, self.children[-1].x2
			
			return xx2 - xx1
	
	def getMergedLeafs(self) :
		#TODO
		pass

if __name__== "__main__" :
	s = SegmentTree()
	s.insert(0, 10, 'region 1')
	s.insert(8, 12, 'region 2')
	s.insert(5, 8, 'region 3')
	s.insert(34, 40, 'region 4')
	s.insert(35, 38, 'region 5')
	s.insert(36, 37, 'region 6')
	print "Tree:"
	print s
	print "removing gaps"
	s.removeGaps()
	#print "flatten"
	#s.flatten()
	print s
	print s.getIndexedLength()
	print s.intersect(16)
	
