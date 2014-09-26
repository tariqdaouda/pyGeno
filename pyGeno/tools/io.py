import sys

def printf(*s) :
	'print + sys.stdout.flush()'
	for e in s[:-1] :
		print e,
	print s[-1]

	sys.stdout.flush()

def enterConfirm_prompt(enterMsg) :
	stopi = False
	while not stopi :
		print "====\n At any time you can quit by entering 'quit'\n===="
		vali = raw_input(enterMsg)
		if vali.lower() == 'quit' :
			vali = None
			stopi = True
		else :
			print "You've entered:\n\t%s" % vali
			valj = confirm_prompt("")
			if valj == 'yes' :
				stopi = True
			if valj == 'quit' :
				vali = None
				stopi = True
				
	return vali

def confirm_prompt(preMsg) :
	while True :
		val = raw_input('%splease confirm ("yes", "no", "quit"): ' % preMsg)
		if val.lower() == 'yes' or val.lower() == 'no' or val.lower() == 'quit':
			return val.lower()
		
