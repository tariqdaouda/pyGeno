import random, copy

def __insertTree(childTree, parentTree):
	"""This a private (You shouldn't have to call it) recursive function that inserts a child tree into a parent tree."""
	if childTree.x1 != None and childTree.x2 != None :
		parentTree.insert(childTree.x1, childTree.x2, childTree.name, childTree.referedObject)

	for c in childTree.children:
		__insertTree(c, parentTree)
		
def __moveTree(newX1, tree):
	"""This a private recursive (You shouldn't have to call it) function that translates tree(and it's children) to a given x1"""
	if tree.x1 != None and tree.x2 != None :
		tree.x1, tree.x2 = newX1, tree.x2-newX1

	for c in childTree.children:
		__moveTree(parentTree, c)
		
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
	
	def __init__(self, x1 = None, x2 = None, name = '', father = None, level = 0, referedObject = None) :
		#print "aaaa"
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
		#print "index", x1, x2, name, self.id
		#print '===', self.name, self.referedObject
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

		#print "NAME", name, self.name, self.id
		rt = None
		insertId = None
		childrenToRemove = []
		for i in range(len(self.children)) :
			if self.children[i].x1 == xx1 and xx2 == self.children[i].x2 :
				return self.children[i]
			
			if self.children[i].x1 <= xx1 and xx2 <= self.children[i].x2 :
				return self.children[i].insert(x1, x2, name, referedObject)
			
			elif xx1 <= self.children[i].x1 and self.children[i].x2 <= xx2 :
				#print "aaaaaaaaaaa", xx1, xx2
				if rt == None :
					rt = SegmentTree(xx1, xx2, name, self, self.level+1, referedObject)
					insertId = i
					
				rt.__addChild(self.children[i])
				#print '===========>', rt
				self.children[i].father = rt
				childrenToRemove.append(self.children[i])
				#childrenToRemove.append(i)
			elif xx1 <= self.children[i].x1 and xx2 <= self.children[i].x2 :
				insertId = i
				break
				
		if rt != None :
			#print 'childrenToRemove', childrenToRemove, len(self.children)
			self.__addChild(rt, insertId)
			for c in childrenToRemove :
				self.children.remove(c)
				#print len(self.children)
		else :
			rt = SegmentTree(xx1, xx2, name, self, self.level+1, referedObject)
			if insertId != None :
				self.__addChild(rt, insertId)
			else :
				self.__addChild(rt)
		
		return rt
	
	def insertTree(self, segTree):
		"""inserts segTree in the right position (regions will rearanged to fit the organisation of self)"""
		__insertTree(self, segTree)
	
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
		
		elif (self.x1 == None and self.x2 == None) : 
			for c in self.children :
				ret.extend(c.intersect(x1, x2))

		return ret

	def bridgeGaps(self) :
		"""returns a tree similar to self but where the gaps between succesive indexed regions have been removed:
		Regions are moved so they adajcent to the one before, the transfornations also affect children"""
		res = SegmentTree()
		
		cc = copy.copy(self.children[0])
		res.insertTree(cc)
		previousC = cc
		for i in range(1, len(self.children)) :
			cc = copy.copy(self.children[i])
			if self.children[i].x1 > previousC.x2:				
				__moveTree(previousC.x2, cc)
			res.insertTree(cc)
			previousC = cc
		return res
		
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
				#print '\n====\n', c
				res.append((c.x1, c.x2)) 
		else :
			if self.x1 != None :
				res = [(self.x1, self.x2)]
			else :
				res = None
		return res

	def flattened(self) :
		r"""Returns a flattend version of the tree. A list of SegmentTrees where all overlapping children have been merged together"""
		#print "Warning: flattened a checker"
		root = SegmentTree()
		#res = []
		if self.x1 == None and self.x2 == None:
			if len(self.children) == 1 :
				#res = [self.children[0]]
				root.insert(self.children[0].x1, self.children[0].x2, '', [self.children[0].referedObject])
			elif len(self.children) > 1 :
				x1 = self.children[0].x1
				x2 = self.children[0].x2
				refObjs = [self.children[0].referedObject]
				for c in self.children[1:]:
					if x2 >= c.x1 :
						#print 'merging', x1, x2, "-", c.x1, c.x2
						x2 = c.x2
						refObjs.append(c.referedObject)
					else :
						# print refObjs
						#res.append(SegmentTree(x1, x2, father = self, referedObject = refObjs))
						root.insert(x1, x2, '', referedObject = refObjs)
						x1 = c.x1
						x2 = c.x2
						refObjs = [c.referedObject]
						
				#res.append(SegmentTree(x1, x2, '', referedObject = refObjs))
				root.insert(x1, x2, '', referedObject = refObjs)
			#else :
			#	res = []
		else :
			#res = [self]
			return self
		#print res
		#return res
		return root
		
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

'''
s = SegmentTree()

s.insert(13, 15)
s.insert(7, 8)
s.insert(9, 10)
s.insert(11, 14)
s.insert(12, 14)
s.insert(0, 12)
s.insert(1, 6)
s.insert(1, 6)
s.insert(1, 6)
s.insert(2, 3)
s.insert(4, 5)

print s

print s.intersect(15, 15) 
'''
