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
	parser.add_option("-l",dest="list",type="string",metavar="FILE",
		help="List of particles (numbered in EMAN format)")
        parser.add_option("-o",dest="out",type="string",metavar="FILE",
                help="Output filename")
        parser.add_option("-b", action="store_true",dest="bad",default=False,
                help="Flag if particle list is for particles to EXCLUDE")
        parser.add_option("-g", action="store_true",dest="good",default=False,
                help="Flag if particle list is for particles to INCLUDE")
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
	list = params['list']	#EMAN numbering scheme
	debug = params['debug']
	fnew = f[:-4]
	bad = params['bad']
	good = params['good']
	if debug is True:
		print 'bad flag is %s' %(bad)
		print 'particle list is %s' %(list)
	fout = params['out']

	o = open(fout,'w')

	f1 = open(f,'r')

	b1 = open(list,'r')
	bad2 = b1.readlines()

	tot = len(f1.readlines())
	tot2 = len(bad2)
	i = 1

	if bad is True:

		if debug is True:
			print 'Bad particle list = %s' %(list)

		while i <= tot:
		
			n = '%s\n' %(i-1)
			if debug is True:
				print 'Checking particle %s' %(n)

			if n in bad2:
				if debug is True:
					print ('---------------------> Particle %s is excluded' %(i-1))
				i = i + 1
				continue

			l = linecache.getline(f,i)
			o.write(l)

			i = i + 1

		o.close()

	if good is True:

		if debug is True:
			print 'Total # particles = %s' %(tot)
		while i <= tot2:

			l = linecache.getline(list,i)
			if debug is True:

				print l

			sel = l.split()
			sel = sel[0]
			if debug is True:
				print sel
			sel = int(sel) + 1

			new = linecache.getline(f,sel)
			o.write(new)
			i = i + 1

		o.close()

if __name__ == "__main__":
     params=setupParserOptions()
     main(params)
