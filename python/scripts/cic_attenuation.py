#!/usr/bin/env python

import argparse
from gemstat.matrix import *

import sys

import scipy as S


def bryan_active_cic2(cic,erk,att):
    #return cic - cic*S.power(1+S.power(erk,-1.0)*S.exp(-att),-1.0)
    return cic*S.power(1.0 + erk*S.exp(att),-1.0)#Identical to above, solved differently.

def hassan_active_cic(cic,erk,cic_att=16.0):
    """Hassan's method of calculating effective CIC concentration."""
    return cic * S.exp(-cic_att * erk)

parser = argparse.ArgumentParser()

parser.add_argument("--C", metavar="C", type=float, help="attenuation parameter.")
parser.add_argument("--target", metavar="TARGET", type=str, help="which factor to attentuate")
parser.add_argument("--attenuator", metavar="SOURCE", type=str, help="name of TF/kinase/input that causes attenuation")

parser.add_argument("--other_matrix", metavar="OTHER", type=str, help="If the attenuator is in a different expression matrix.")

parser.add_argument("--bryan",action="store_true", help="Use bryan's formula for CIC attenuation, instead of hassan's.")

parser.add_argument("--prescale_attenuator",metavar="P",type=float,default=1.0,help="Scale the attenuator before applying attenuation.")

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

attenuator *= args.prescale_attenuator

after_attenuation = None
if args.bryan:
	after_attenuation = bryan_active_cic2(in_matrix.storage[target_row_id,:], attenuator, args.C)
else:
	after_attenuation = hassan_active_cic(in_matrix.storage[target_row_id,:], attenuator, args.C)


in_matrix.storage[target_row_id,:] = after_attenuation

in_matrix.write(args.OUTFILE,format="%.5f")
