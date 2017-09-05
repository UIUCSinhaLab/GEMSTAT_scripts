#!/usr/bin/env python

import argparse
from gemstat.matrix import *

import sys

import scipy as S

parser = argparse.ArgumentParser()

parser.add_argument("--S", metavar="S", type=float, default=1.0, help="scale parameter.")
parser.add_argument("--target", metavar="TARGET", type=str, default=None, help="which factor to scale (empty scales all)")

parser.add_argument("INFILE", metavar="IN_FILE", type=str)
parser.add_argument("OUTFILE", metavar="OUT_FILE", type=str)

args, other = parser.parse_known_args()


in_matrix = GEMSTAT_Matrix.load(args.INFILE)

if None == args.target: #scale everything
	in_matrix.storage *= args.S
else:#scale a specific entry only
	target_row_id = in_matrix.names_to_rows[args.target]
	in_matrix.storage[target_row_id,:] *= args.S

in_matrix.write(args.OUTFILE,format="%.5f")
