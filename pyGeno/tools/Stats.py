import numpy as np

def kullback_leibler(p, q) :
	"""Discrete Kullback-Leibler divergence D(P||Q)"""
	p = np.asarray(p, dtype=np.float)
	q = np.asarray(q, dtype=np.float)

	if p.shape != q.shape :
		raise ValueError("p and q must be of the same dimensions")
	
	return np.sum(np.where(p > 0, np.log(p / q) * p, 0))

def squaredError_log10(p, q) :
	p = np.asarray(p, dtype=np.float)
	q = np.asarray(q, dtype=np.float)
	
	if p.shape != q.shape :
		raise ValueError("p and q must be of the same dimensions")
		
	return np.log10(sum((p-q)**2)) - np.log(len(p))
	
def fisherExactTest(table) :
	"""Fisher's exact test
	----------
	table: contengency table
	"""
	raise NotImplementedError
