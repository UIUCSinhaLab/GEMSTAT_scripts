#!/usr/bin/env python
from __future__ import print_function

import scipy as S
import sys

import re

sys.stderr.write("NEVER EVER EVER run this on a server, it executes arbitrary code provided in the template file.\n")

def ret_exec(in_string):
	retval = None
	exec("retval = %s" % in_string)
	return retval

class mydict(dict):
    def __getitem__(self,key):
        if isinstance(key,str) and key.find(":") > -1:
            f_name,the_args = key.split(":")
            retval = None
            exec("retval = %s(%s)" %( f_name, the_args))
            return retval
        else:
            return super(mydict, self).__getitem__(key)


def uniform(low,high):
	return S.random.uniform(low,high)

def log_uniform(low,high):
	return S.exp(uniform(S.log(low),S.log(high)))

def const(val):
	return val

def log(low,high):
	return S.log(uniform(S.exp(low),S.exp(high)))

def normal(mu,std):
	return S.random.normal(mu,std)

#K_range = S.log(S.array([0.01, 10000]))
#a_act_range = S.log(S.array([1, 10]))
#a_rep_range = S.log(S.array([0.00001, 1]))
#coop_range = S.array([1, 100])
#q_btm_range = S.array([0.001, 0.01])

import pystache

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("INFILE", metavar="IN_FILE", type=str)
#parser.add_argument("OUTFILE", metavar="OUT_FILE", type=str)
args, other = parser.parse_known_args()

foo = mydict()

thetemplate = open(args.INFILE).read()

my_regex = re.compile("{{(.*?)}}")

things = my_regex.findall(thetemplate)
replaced_things = [foo[one_thing] for one_thing in things]

final_template = my_regex.sub("%f",thetemplate)
print(final_template.strip() % tuple(replaced_things))
