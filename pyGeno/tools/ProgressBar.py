import sys, time, cPickle

class ProgressBar :
	"""A very simple unthreaded progress bar. This progress bar also logs stats in .logs.
	Usage example::

		p = ProgressBar(nbEpochs = -1)
			for i in range(200000) :
				p.update(label = 'value of i %d' % i)
		p.close()
	
	If you don't know the maximum number of epochs you can enter nbEpochs < 1
	"""
	
	def __init__(self, nbEpochs = -1, width = 25, label = "progress", minRefeshTime = 1) :
		self.width = width
		self.currEpoch = 0
		self.nbEpochs = float(nbEpochs)
		self.bar = ''

		self.label = label
		self.wheel = ["-", "\\", "|", "/"]
		self.startTime = time.time()
		self.lastPrintTime = -1
		self.minRefeshTime = minRefeshTime
		
		self.runtime = -1
		self.runtime_hr = -1
		self.avg = -1
		self.remtime = -1
		self.remtime_hr = -1
		self.currTime = time.time()
		self.lastEpochDuration = -1 
		
		self.bars = []
		self.miniSnake = '~-~-~-?:>' 
		self.logs = {'epochDuration' : [], 'avg' : [], 'runtime' : [], 'remtime' : []}
		
	def formatTime(self, val) :
		if val < 60 :
			return '%.3fsc' % val
		elif val < 3600 :
			return '%.3fmin' % (val/60)
		else :
			return '%dh %dmin' % (int(val)/3600, int(val/60)%60)

	def _update(self) :
		tim = time.time()
		if self.nbEpochs > 1 :
			if self.currTime > 0 :
				self.lastEpochDuration = tim - self.currTime
			self.currTime = tim
			self.runtime = tim - self.startTime
			self.runtime_hr = self.formatTime(self.runtime)
			self.avg = self.runtime/self.currEpoch
			
			self.remtime = self.avg * (self.nbEpochs-self.currEpoch)
			self.remtime_hr = self.formatTime(self.remtime)
	
	def log(self) :
		"""logs stats about the progression, without printing anything on screen"""
		
		self.logs['epochDuration'].append(self.lastEpochDuration)
		self.logs['avg'].append(self.avg)
		self.logs['runtime'].append(self.runtime)
		self.logs['remtime'].append(self.remtime)
	
	def saveLogs(self, filename) :
		"""dumps logs into a nice pickle"""
		f = open(filename, 'wb')
		cPickle.dump(self.logs, f)
		f.close()

	def update(self, label = '', forceRefresh = False, log = False) :
		"""the function to be called at each iteration. Setting log = True is the same as calling log() just after update()"""
		self.currEpoch += 1
		tim = time.time()
		if (tim - self.lastPrintTime > self.minRefeshTime) or forceRefresh :
			self._update()
			
			wheelState = self.wheel[self.currEpoch%len(self.wheel)]
			
			if label == '' :
				slabel = self.label
			else :
				slabel = label
			
			if self.nbEpochs > 1 :
				ratio = self.currEpoch/self.nbEpochs
				snakeLen = int(self.width*ratio)
				voidLen = int(self.width - (self.width*ratio))

				if snakeLen + voidLen < self.width :
					snakeLen = self.width - voidLen
				
				self.bar = "%s %s[%s:>%s] %.2f%% (%d/%d) runtime: %s, remaing: %s, avg: %s" %(wheelState, slabel, "~-" * snakeLen, "  " * voidLen, ratio*100, self.currEpoch, self.nbEpochs, self.runtime_hr, self.remtime_hr, self.formatTime(self.avg))
				if self.currEpoch == self.nbEpochs :
					self.close()
				
			else :
				w = self.width - len(self.miniSnake)
				v = self.currEpoch%(w+1)
				snake = "%s%s%s" %("  " * (v), self.miniSnake, "  " * (w-v))
				self.bar = "%s %s[%s] %s%% (%d/%s) runtime: %s, remaining: %s, avg: %s" %(wheelState, slabel, snake, '?', self.currEpoch, '?', self.runtime_hr, '?', self.formatTime(self.avg))
			
			sys.stdout.write("\b" * (len(self.bar)+1))
			sys.stdout.write(" " * (len(self.bar)+1))
			sys.stdout.write("\b" * (len(self.bar)+1))
			sys.stdout.write(self.bar)
			sys.stdout.flush()
			self.lastPrintTime = time.time()
			
			if log :
				self.log()

	def close(self) :
		"""Closes the bar so your next print will be on another line"""
		self.update(forceRefresh = True)
		print '\n'
		
if __name__ == "__main__" :
	p = ProgressBar(nbEpochs = 100000000000)
	for i in xrange(100000000000) :
		p.update()
		#time.sleep(3)
	p.close()
	
