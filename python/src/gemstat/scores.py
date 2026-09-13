import numpy as _NP
#import scipy as _S

def sse_old(gt,prediction):
	"""Weighting not yet implemented.
	"""

	return _NP.real(_NP.power(gt - prediction,2.0).sum())

def sse(gt,predictions):
	gt = _NP.array(gt,ndmin=2)
	predictions = _NP.array(predictions,ndmin=2)
	N,M = predictions.shape
	errors = _NP.power(_NP.maximum(0.0,predictions) - _NP.maximum(0.0,gt), 2.0).mean(1)
	return _NP.absolute(_NP.real(errors))

def wPGP(gt,prediction):
	"""
	Not yet weighted.
	"""
	r_max = gt.max()

	prediction = _NP.minimum(prediction,1.0)

	reward = (gt*_NP.minimum(gt,prediction)).sum()/_NP.power(gt,2.0).sum()

	penalty_num = (
		(r_max - gt)
		*(prediction - gt)
		*_NP.array(prediction > gt,dtype=_NP.float_)
		).sum()

	penalty_denom = _NP.power(r_max - gt,2.0).sum()
	penalty = penalty_num / penalty_denom

	wpgp_score = 0.5 + 0.5*(reward-penalty)
	assert wpgp_score <= 1.0 and wpgp_score >= 0.0, "PGP score not in acceptable range"
	return wpgp_score
