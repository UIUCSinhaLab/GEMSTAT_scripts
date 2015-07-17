def sse(gt,prediction):
	"""Weighting not yet implemented.
	"""
	
	return S.power(gt - prediction).mean()

def wPGP(gt_prediction):
	"""
	Not yet weighted.
	"""
	r_max = gt.max()
	
	prediction = S.minimum(prediction,1.0)
	
	reward = (gt*S.minimum(gt,prediction)).sum()/S.power(gt,2.0).sum()
	
    	penalty_num = (
		(r_max - gt)
		*(prediction - gt)
		*S.array(prediction > gt,dtype=S.float_)
		).sum()
	
	penalty_denom = S.power(r_max - gt,2.0).sum()
	penalty = penalty_num / penalty_denom
	
	wpgp_score = 0.5 + 0.5*(reward-penalty)
	assert wpgp_score <= 1.0 and wpgp_score >= 0.0, "PGP score not in acceptable range"
	return wpgp_score
