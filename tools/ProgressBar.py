import sys, time

class ProgressBar :
	def __init__(self, nbEpochs, width = 50, name = "progress") :
		self.width = width
		self.currEpoch = 0
		self.nbEpochs = float(nbEpochs)
		self.bar = ''
		
		self.name = name
		self.wheel = ["-", "\\", "|", "/"]
		
	def update(self) :
		self.currEpoch += 1
		ratio = self.currEpoch/self.nbEpochs
		snakeLen = int(self.width*ratio)
		voidLen = int(self.width - (self.width*ratio))
		
		wheelState = self.wheel[self.currEpoch%len(self.wheel)]
		
		if snakeLen + voidLen < self.width :
			snakeLen = self.width - voidLen
	
		self.bar = "%s %s[%s:>%s] %.1f%%" %(wheelState, self.name, "~-" * snakeLen, "  " * voidLen, ratio*100)
		sys.stdout.write(self.bar)
		sys.stdout.flush()
		sys.stdout.write("\b" * (len(self.bar)+1))

if __name__ == '__main__' :
	p = ProgressBar(nbEpochs = 200)
	for i in range(200) :
		p.update()
		time.sleep(0.1)
