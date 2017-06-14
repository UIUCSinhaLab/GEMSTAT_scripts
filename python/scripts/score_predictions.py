#!/usr/bin/env python

import argparse
from gemstat.matrix import *
import gemstat.scores as scoring

import sys

import scipy as S

parser = argparse.ArgumentParser()

parser.add_argument("-n",default=False,action="store_true",help="Don't print the newline after the score. Sometimes useful.")
parser.add_argument("--score",metavar="SCORE",type=str,help="One of SSE, PGP")

parser.add_argument("--GT_matrix", metavar="GT",default=None, type=str, help="Use this to substitute the values out of this matrix for the ground truth, instead of using that out of the file to be scored.")
parser.add_argument("INFILE", metavar="IN_FILE", type=str)

args, other = parser.parse_known_args()


in_matrix = GEMSTAT_Matrix.load(args.INFILE)
N, M = in_matrix.storage.shape

other_mat = None
N_other = None
if args.GT_matrix:
	other_mat = GEMSTAT_Matrix.load(args.GT_matrix)
	N_other, _ = other_mat.storage.shape
	print "loaded other matrix"

predictions = None
ground_truth = None

if N_other == None:
	ground_truth, predictions = in_matrix.separate_output()
elif N_other == N/2:
	ground_truth = other_mat
	_,predictions = in_matrix.separate_output()
elif N_other == N:
	ground_truth = other_mat
	predictions = in_matrix
else:
	assert False, "This should not happen"
#End of ifs

the_score = None

if args.score == "SSE":
	the_score = scoring.sse(ground_truth.storage,predictions.storage).mean()
elif args.score == "PGP":
	the_score = scoring.wPGP(ground_truth.storage,predictions.storage)
else:
	assert False, "you didn't provide a valid scoring function"


sys.stdout.write("%.10f" % the_score)
if not args.n:
	sys.stdout.write("\n")
sys.stdout.flush()
