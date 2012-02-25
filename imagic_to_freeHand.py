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
	parser.set_usage("%prog -f <stack> -s <select> --box=<float>")
	parser.add_option("-f",dest="stack",type="string",metavar="FILE",
		help="raw, IMAGIC particle stack from tilted data (black particles)")
        parser.add_option("-s",dest="sel",type="string",metavar="FILE",
                help="Select file generated b 'EMAN_to_FREALIGN_freeHand.py'")
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

def getIMAGICPath():
        ### get the imagicroot directory
        impath = subprocess.Popen("env | grep IMAGIC_ROOT", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        imagicpath = impath.replace("IMAGIC_ROOT=","")
        if imagicpath != '/opt/qb3/imagic-070813':
                        print "imagic/070813 was not found, make sure it is in your path"
                        sys.exit()

def getImagicVersion(imagicroot):
        ### get IMAGIC version from the "version_######S" file in
        ### the imagicroot directory, return as an int
        versionstr=glob.glob(os.path.join(imagicroot,"version_*"))
        if versionstr:
                v = re.search('\d\d\d\d\d\d',versionstr[0]).group(0)
                return int(v)
        else:
                print "Could not get version number from imagic root directory"
		sys.exit()

def main(params):
	
	stack = params['stack']
	base = stack[:-4]
	box = str(params['box'])
	sel = params['sel']

	f = open(sel,'r')
	tot = len(f.readlines())
	f.close()

	cmd="e2proc2d.py --list=%s %s.img %s_model1.img" %(sel,base,base)
	subprocess.Popen(cmd,shell=True).wait()

	cmd="~michael/BATCHLIB/freeHand/testim.b %s %s" %(box[:-2],tot)
	subprocess.Popen(cmd,shell=True).wait()

	cmd="rm %s_model1.hed" %(base)
	subprocess.Popen(cmd,shell=True).wait()

	cmd="cp test.hed %s_model1.hed" %(base)
	subprocess.Popen(cmd,shell=True).wait()

	cmd="~michael/BATCHLIB/freeHand/imagic_to_3Dmrc.b %s_model1.img" %(base)
	subprocess.Popen(cmd,shell=True).wait()

if __name__ == "__main__":
     imagicroot = getIMAGICPath()
     params=setupParserOptions()
     main(params)

