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
	parser.set_usage("%prog -p <parameter> -b <exclude>")
	parser.add_option("-p",dest="param",type="string",metavar="FILE",
		help="Parameter file with per-particle CTF information from ctftilt.py")
	parser.add_option("-b",dest="bad",type="string",metavar="FILE",
		help="List of particles to exclude (numbered in EMAN format)")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
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

#Exclude parameters from list

#Inputs 
def main(params):
	f = params['param']
	bad = params['bad']	#EMAN numbering scheme
	debug = params['debug']
	fnew = f[:-4]

	fout = '%s_sel.par' %(fnew)

	o = open(fout,'w')

	f1 = open(f,'r')

	b1 = open(bad,'r')
	bad = b1.readlines()

	tot = len(f1.readlines())

	i = 1

	while i <= tot:
	
		n = '%s\n' %(i-1)
	
		if n in bad:
			if debug is True:
				print ('Particle %s is excluded' %(i))
			i = i + 1
			continue

		l = linecache.getline(f,i)
		line = l.split()
		o.write('%s		%s		%s		%s\n' %(line[0],line[1],line[2],line[3]))

		i = i + 1

	o.close()

if __name__ == "__main__":
     params=setupParserOptions()
     main(params)
