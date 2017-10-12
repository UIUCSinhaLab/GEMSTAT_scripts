#!/usr/bin/env python

import argparse
from gemstat.matrix import *

import sys

import scipy as S

parser = argparse.ArgumentParser()
parser.add_argument('--s0', metavar="S0", type=float, nargs="+", help="Zeroeth coefficient for stddev polynomial.")
parser.add_argument("--s1", metavar="S1", type=float, nargs="+", help="First coefficient for stddev polynomial. (Stddev will be S0 + S1*(mu_bin) for each bin")

parser.add_argument("--baseshift", metavar="BASESHIFT", type=float, nargs="+", help="shift/scale curves so that the 0.0 line shifts to this value, 1.0 line remains the same.")

parser.add_argument("--C", metavar="C", type=int, help="Number of expanded curves to include")
parser.add_argument("INFILE", metavar="IN_FILE", type=str)
parser.add_argument("OUTFILE", metavar="OUT_FILE", type=str)

args, other = parser.parse_known_args()

in_matrix = GEMSTAT_Matrix.load(args.INFILE)

N, M = in_matrix.storage.shape
MAXES = in_matrix.storage.max(1)

S0 = S.array(args.s0 if len(args.s0) == N else S.array([args.s0[0] for i in range(N)]))
S1 = S.array(args.s1 if len(args.s1) == N else S.array([args.s1[0] for i in range(N)]))

if None != args.baseshift and len(args.baseshift) > 0:
	if len(args.baseshift) == 1:
		args.baseshift = args.baseshift * N
	for i in range(N):
		in_matrix.storage[i,:] = in_matrix.storage[i,:] * (1.0 - args.baseshift[i]/MAXES[i]) + args.baseshift[i]

		
	 



def basic_variance(curve, s0, s1):
	return s0 + s1*curve


CHOSEN_VARIANCE=basic_variance

copies = list()
for copy_no in range(args.C):
	fuzzed_rows = list()
	for i in range(N):
		row = in_matrix.storage[i]
		fuzzed = None #I know python scoping doesn't require this, but it is more understandable.
		if S0[i] == 0.0 and S1[i] == 0.0:
			fuzzed = row
		else:
			fuzzed = S.random.normal(loc=row, scale=CHOSEN_VARIANCE(row, S0[i], S1[i]))
			fuzzed = S.maximum(0.0, fuzzed) # noized input cannot be less than 0.0
			#fuzzed = S.minimum(1.0, fuzzed) # Does it make sense to constrain the top end? For now, no. Placeholder code.
		fuzzed_rows.append(fuzzed)
	copies.append(S.array(fuzzed_rows))
in_matrix.storage = S.hstack([in_matrix.storage] + copies)

in_matrix.write(args.OUTFILE)
