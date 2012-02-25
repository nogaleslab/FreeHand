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
        parser.set_usage("%prog -i <input> -o <output> --pix=<float>")
        parser.add_option("-i",dest="input",type="string",metavar="FILE",
                help="FREALIGN parameter input file")
        parser.add_option("-o",dest="output",type="string",metavar="FILE",
                help="Reformatted output parameter filename")
        parser.add_option("--pix",dest="pix",type="float", metavar="INT",
                help="Pixel size")
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

	f = params['input']
	o = params['output']

	f1 = open(f,'r')
	pix = float(params['pix'])
	o1 = open(o,'wa')

	o1.write("C Frealign format parameter file created from Search_fspace parameter file\n")
	o1.write("C\n")
	o1.write("C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMax\n")

	for line in f1:

		l = line.split()
		if l[0] == 'C':
			continue
		o1.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.0f%6d%9.1f%9.1f%8.2f%7.2f\n" %(float(l[0]),float(l[1]),float(l[2]),float(l[3]),float(l[4])/pix,float(l[5])/pix,float(l[6]),float(l[7]),float(l[8]),float(l[9]),float(l[10]),float(l[13])))

	o1.write("C\n")

if __name__ == "__main__":
     params=setupParserOptions()
     main(params)
