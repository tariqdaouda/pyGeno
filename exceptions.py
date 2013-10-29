class ConfigurationError(Exception) :
	def __init__(self, msg) :
		self.msg = msg
	def __str__(self) :
		return self.msg + '\n'
		
class ChromosomeError(Exception) :
	def __init__(self, msg, number) :
		self.msg = msg
		self.number = number
	def __str__(self) :
		return self.msg + '\n'

class GenomeError(Exception) :
	def __init__(self, msg, path) :
		self.msg = msg
		self.path = path
	def __str__(self) :
		return self.msg + '\n'

class GeneNotFound(Exception):
	def __init__(self, chromosome, geneSymbol, message = ''):
		self.message = message
		self.symbol = geneSymbol
		self.chromosome = chromosome
		
	def __str__(self):
		return """
		Description : %s
		gene_symbol : %s
		chromosome : %s\n"""%(self.message, self.symbol, self.chromosome)

class SNPError(Exception) :
	def __init__(self, message, line):
		self.message = message
		self.line = line
	
	def __str__(self):
		return "Invalid SNP, %s" % self.message

class DebugException(Exception) :
	"A dummy multi purpose exception"
	def __init__(self, msg, moreData) :
		self.msg = msg
		self.moreData = moreData
	def __str__(self) :
		return self.msg + '\n'


class RequestError(Exception):
	def __init__(self, message = ''):
		self.message = message
		
	def __str__(self):
		return """Request Error: %s"""%(self.message)
