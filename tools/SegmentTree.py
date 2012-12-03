import random

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
			rt = SegmentTree(xx1, xx2, name, self, referedObject)
			if insertId != None :
				self.__addChild(rt, insertId)
			else :
				self.__addChild(rt)
		
		return rt
		
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

	def getX1(self) :
		if self.x1 != None :
			return self.x1
		return self.children[0].x1

	def getX2(self) :
		if self.x2 != None :
			return self.x2
		return self.children[-1].x2
	
	def getIndexedRegionsLength(self) :
		if self.x1 != None :
			return self.x2 - self.x1
		else :
			if len(self.children) == 0 :
				return 0
			else :
				l = 0
				for c in self.children :
					l += c.x2 - c.x1
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

	def getMergedFirstLevel(self) :
		'same as getFirstLevel, but overlapping children are merged'
		#print "Warning: getMergedFirstLevel a checker"
		res = []
		if len(self.children) == 1 :
			res.append((self.children[0].x1, self.children[0].x2))
		elif len(self.children) > 1 :
			x1 = self.children[0].x1
			x2 = self.children[0].x2
			for c in self.children[1:]:
				if x2 > c.x1 :
					x2 = c.x2
				else :
					res.append((x1, x2))
					x1 = c.x1
					x2 = c.x2
			res.append((x1, x2))
		else :
			if self.x1 != None :
				res = [(self.x1, self.x2)]
			else :
				res = None
		return res
		
	def __str__(self) :
		strRes = repr(self)
		
		offset = ''
		for i in range(self.level+1) :
			offset += '\t'
			
		for c in self.children :
			strRes += '\n'+offset+'-->'+str(c)
		
		return strRes
	
	def __repr__(self) :
		if self.x1 == None :
			if len(self.children) > 0 :
				return "Root : %d-%d, name : %s, id : %d" %(self.children[0].x1, self.children[-1].x2, self.name, self.id)
			else :
				return "Root : EMPTY , name : %s, id : %d" %(self.name, self.id)
		else :
			return "Segment : %d-%d, name : %s, id : %d, father id : %d" %(self.x1, self.x2, self.name, self.id, self.father.id)
			
		
	def __len__(self) :
		r"""if root returns whole region of interest length. if not it is identical to getIndexedRegionsLength()"""
		
		if self.x1 != None :
			return self.getIndexedRegionsLength()
		else :
			if len(self.children) == 0 :
				return 0
			else :
				xx1, xx2 = self.children[0].x1, self.children[-1].x2
			
			return xx2 - xx1
	
	def getEffectiveLength(self) :
		r"returns the sum of the total length of the leafs without overlap"
		#print "getEffectiveLength a cheker"
		if self.x1 == None and len(self.children) == 0 :
			return 0
			
		if len(self.children) == 0 :
			return self.x2 - self.x1
		else :
			if len(self.children) == 1 :
				return self.children[0].x2 - self.children[0].x1
			else :
				l = 0
				for i in range(0, len(self.children)-1) :
					l += self.children[i].getEffectiveLength() - max(self.children[i].x2 - self.children[i+1].x1, 0)		
				return l+self.children[-1].getEffectiveLength()
	
	def getMergedLeafs(self) :
		#TODO
		pass
	
	def getIndexedRegionsLength(self) :
		if self.x1 != None :
			return self.x2 - self.x1
		else :
			if len(self.children) == 0 :
				return 0
			else :
				return self.children[-1].x2 - self.children[0].x1

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
