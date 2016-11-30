import scipy as _S

def sse(gt,prediction):
	"""Weighting not yet implemented.
	"""
	
	return _S.real(_S.power(gt - prediction,2.0).mean())

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
