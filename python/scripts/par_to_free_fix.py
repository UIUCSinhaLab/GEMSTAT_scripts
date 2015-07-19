#!/usr/bin/env python
from __future__ import print_function

import sys
import argparse

def read_factor_info(file):
	lines = [line.strip().split() for line in file]
	return dict([(l[0],(int(l[1]),int(l[2]))) for l in lines])


from gemstat.model import ParFile, parfile_to_vector

parser = argparse.ArgumentParser()
parser.add_argument("-f","--factor_info",dest="factor_info",type=str)
parser.add_argument("INFILE", metavar="IN_FILE", type=str)
parser.add_argument("OUTFILE", default="-", metavar="OUT_FILE", type=str)
args, other = parser.parse_known_args()

#IN_FILE = open(args.IN_FILE)
OUT_FILE = sys.stdout if args.OUTFILE == "-" else open(args.OUTFILE)

factor_info = None
if args.factor_info:
	fifile = open(args.factor_info)
	factor_info = read_factor_info(fifile)
	fifile.close()


dict_rep = ParFile.read_old_par(args.INFILE)
vector = parfile_to_vector(dict_rep,factor_roles=factor_info,include_pis=True)

for i in vector:
	print(int(i),file=OUT_FILE)
