#!/usr/bin/env python

import optparse
from sys import *
import os,sys,re
from optparse import OptionParser
import glob
import subprocess
from os import system
import linecache

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -f <stack> --total= <float> --box=<float>")
	parser.add_option("-f",dest="stack",type="string",metavar="FILE",
		help="raw, IMAGIC particle stack from tilted data (black particles)")
        parser.add_option("--total",dest="tot",type="float", metavar="FLOAT",
                help="Number of particles")
	parser.add_option("--box",dest="box",type="float", metavar="FLOAT",
		help="Box size")
	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit()
	params={}
	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params

def main(params):
	
	stack = params['stack']
	base = stack[:-4]
	box = str(params['box'])
	tot = str(params['tot'])
	
	cmd="~michael/BATCHLIB/freeHand/testim.b %s %s" %(box[:-2],tot[:-2])
	subprocess.Popen(cmd,shell=True).wait()

	cmd="rm %s.hed" %(base)
	subprocess.Popen(cmd,shell=True).wait()

	cmd="cp test.hed %s.hed" %(base)
	subprocess.Popen(cmd,shell=True).wait()

	cmd="~michael/BATCHLIB/freeHand/imagic_to_3Dmrc.b %s.img" %(base)
	subprocess.Popen(cmd,shell=True).wait()

if __name__ == "__main__":
     params=setupParserOptions()
     main(params)

