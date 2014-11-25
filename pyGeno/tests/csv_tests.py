import unittest
from pyGeno.tools.parsers.CSVTools import *

class CSVTests(unittest.TestCase):
		
	def setUp(self):
		pass

	def tearDown(self):
		pass

	def test_createParse(self) :
		testVals = ["test", "test2"]
		c = CSVFile(legend = ["col1", "col2"], separator = "\t")
		l = c.newLine()
		l["col1"] = testVals[0]
		l = c.newLine()
		l["col1"] = testVals[1]
		c.save("test.csv")
		
		c2 = CSVFile()
		c2.parse("test.csv", separator = "\t")
		i = 0
		for l in c2 :
			self.assertEqual(l["col1"], testVals[i])
			i += 1

def runTests() :
	unittest.main()

if __name__ == "__main__" :
	runTests()
