#!/usr/bin/env python

import argparse
from gemstat.matrix import *

import sys

import scipy as S

parser = argparse.ArgumentParser()

parser.add_argument("--C", metavar="C", type=float, help="attenuation parameter.")
parser.add_argument("--target", metavar="TARGET", type=str, help="which factor to attentuate")
parser.add_argument("--attenuator", metavar="SOURCE", type=str, help="name of TF that causes attenuation")

parser.add_argument("--other_matrix", metavar="OTHER", type=str, help="If the attenuator is in a different expression matrix.")

parser.add_argument("INFILE", metavar="IN_FILE", type=str)
parser.add_argument("OUTFILE", metavar="OUT_FILE", type=str)

args, other = parser.parse_known_args()


in_matrix = GEMSTAT_Matrix.load(args.INFILE)
N, M = in_matrix.storage.shape

other_mat = None
if args.other_matrix:
	other_mat = GEMSTAT_Matrix.load(args.other_matrix)
	print "loaded other matrix"

target_row_id = -1
try:
	target_row_id = in_matrix.names_to_rows[args.target]
except:
	pass

if target_row_id < 0:
	print "Could not find the target in the matrix provided."
	sys.exit(1)

attenuator = None

if other_mat:
	attenuator = other_mat.storage[other_mat.names_to_rows[args.attenuator],:]
else:
	attenuator = in_matrix.storage[in_matrix.names_to_rows[args.attenuator],:]

attenuation = S.exp( - args.C * attenuator )

in_matrix.storage[target_row_id,:] *= attenuation

in_matrix.write(args.OUTFILE,format="%.5e")
