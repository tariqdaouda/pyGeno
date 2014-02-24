import sys, time

class ProgressBar :
	"A very simple unthreaded progress bar, see ProgressBar  __name__ == '__main__' for an ex of utilisation"
	def __init__(self, nbEpochs, width = 25, label = "progress", minRefeshTime = 0.1) :
		self.width = width
		self.currEpoch = 0
		self.nbEpochs = float(nbEpochs)
		self.bar = ''

		self.label = label
		self.wheel = ["-", "\\", "|", "/"]
		self.startTime = time.time()
		self.lastPrintTime = 0
		self.minRefeshTime = minRefeshTime

	def formatTime(self, val) :
		if val < 60 :
			return '%.1fsc' % val
		elif val < 3600 :
			return '%.1fmin' % (val/60)
		else :
			return '%dh %dmin' % (int(val)/3600, int(val/60)%60)

	def update(self, label = '') :
		self.currEpoch += 1
		if time.time() - self.lastPrintTime > self.minRefeshTime :
			ratio = self.currEpoch/self.nbEpochs
			snakeLen = int(self.width*ratio)
			voidLen = int(self.width - (self.width*ratio))

			wheelState = self.wheel[self.currEpoch%len(self.wheel)]

			if snakeLen + voidLen < self.width :
				snakeLen = self.width - voidLen

			if label == '' :
				slabel = self.label
			else :
				slabel = label

			elTime = time.time() - self.startTime
			runtime = self.formatTime(elTime)
			remtime = self.formatTime(elTime/self.currEpoch * (self.nbEpochs-self.currEpoch))

			self.bar = "%s %s[%s:>%s] %.1f%% (%d/%d) runtime: %s, remaing: %s" %(wheelState, slabel, "~-" * snakeLen, "  " * voidLen, ratio*100, self.currEpoch, self.nbEpochs, runtime, remtime)
			sys.stdout.write("\b" * (len(self.bar)+1))
			sys.stdout.write(" " * (len(self.bar)+1))
			sys.stdout.write("\b" * (len(self.bar)+1))
			sys.stdout.write(self.bar)
			sys.stdout.flush()
			self.lastPrintTime = time.time()

		if self.currEpoch == self.nbEpochs :
			print

if __name__ == '__main__' :
	p = ProgressBar(nbEpochs = 200)
	for i in range(200) :
		p.update(label = 'value of i %d' % i)
