import scipy as _S

def sse_old(gt,prediction):
	"""Weighting not yet implemented.
	"""
	
	return _S.real(_S.power(gt - prediction,2.0).sum())

def sse(gt,predictions):
    predictions = _S.array(predictions,ndmin=2)
    N,M = predictions.shape
    
    errors = [ _S.power(_S.maximum(0.0,row)-_S.maximum(0.0,gt),2.0).mean(1) for row in predictions]#we expect a vector, but if we get a matrix, this mean is more deffensive. (will cause an error later)
    return _S.absolute(_S.array(errors).ravel())

def wPGP(gt,prediction):
	"""
	Not yet weighted.
	"""
	r_max = gt.max()
	
	prediction = _S.minimum(prediction,1.0)
	
	reward = (gt*_S.minimum(gt,prediction)).sum()/_S.power(gt,2.0).sum()
	
    	penalty_num = (
		(r_max - gt)
		*(prediction - gt)
		*_S.array(prediction > gt,dtype=_S.float_)
		).sum()
	
	penalty_denom = _S.power(r_max - gt,2.0).sum()
	penalty = penalty_num / penalty_denom
	
	wpgp_score = 0.5 + 0.5*(reward-penalty)
	assert wpgp_score <= 1.0 and wpgp_score >= 0.0, "PGP score not in acceptable range"
	return wpgp_score
