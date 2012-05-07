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
        parser.set_usage("%prog -s <stack> -c <ctf info> -p <param>")
        parser.add_option("-s",dest="stack",type="string",metavar="FILE",
                help="Imagic particle stack")
        parser.add_option("-c",dest="ctf",type="string", metavar="FILE",
                help="Per-particle CTF info")
	parser.add_option("-p",dest="param",type="string", metavar="FILE",
                help="Free-hand parameter file (free_param.par)")
        parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 4:
                parser.print_help()
                sys.exit()
        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params

def flip(stack,ctf,boxSize,numParts,Cs,pixSize,amp):
	spifile = "currentSpiderScript.spi"
        if os.path.isfile(spifile):
                os.remove(spifile)
        spi=open(spifile,'w')
	spicmd='x71 = %s\n' %(boxSize)
	spicmd+='x70 = %s\n' %(numParts)
	spicmd+=';CS (MM)\n'
	spicmd+='X50=%s\n' %(Cs)
	spicmd+=';LAMBDA(A)\n'
	spicmd+='X51=0.0336\n'                     
	spicmd+=';WINDOW SIZE(PIX)\n'
	spicmd+='X52=%s\n' %(boxSize)              
	spicmd+=';MAXIMUM SPATIAL FREQUENCY [1/A]\n'
	spicmd+='X53=%s\n' %(str(1/(2*pixSize)))
	spicmd+=';SOURCE SIZE[1/A]\n'
	spicmd+='X54=0.0047\n'
	spicmd+=';DEFOCUS SPREAD[A]\n'
	spicmd+='X55=100\n'
	spicmd+=';ASTIGMATISM[A]\n'
	spicmd+=';X56=0.0\n'
	spicmd+=';AZIMUTH[DEG]\n'
	spicmd+=';X57=0.0\n'
	spicmd+=';AMPLITUDE CONTRAST RATIO [0-1]\n'
	spicmd+='X58=%s\n' %(amp)
	spicmd+='; GAUSSIAN ENVELOPE HALFWIDTH[1/A]\n'
	spicmd+='x59=0.15\n'
	spicmd+=';Sign (+1 or -1)\n'
	spicmd+='X60=-1                  ;-1 to keep the same contrast as input stack\n'
	spicmd+='UD N x70\n'
	spicmd+='%s\n' %(ctf)
	spicmd+='do lb1 x10=1,x70\n'
        spicmd+='UD IC x10,x23,x92,x93\n'
        spicmd+='%s\n' %(ctf)
        spicmd+='x92=x92/2\n'
        spicmd+='TF CT\n'
        spicmd+='%s@{*****x10}\n' %(stack)
        spicmd+='X50             ; CS[mm]\n'
        spicmd+='X23,X51         ; defocus, lambda\n'
        spicmd+='X52             ; dimensions of output array (box size)\n'
        spicmd+='X53             ; max spatial freq\n'
        spicmd+='X54,X55         ; source size, defocus spread\n'
        spicmd+='x92,x93         ; astigmatism correction\n'
        spicmd+='X58,x59         ; amp contrast\n'
	spicmd+='X60             ; sign\n'
        spicmd+='lb1\n'
	spicmd+='UD ICE\n'
	spicmd+='%s\n' %(ctf)
        runSpider(spicmd)

def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            return match.string

def mult(ctfStack,stack,outStack,numParts):
        spifile = "currentSpiderScript.spi"
        if os.path.isfile(spifile):
                os.remove(spifile)
        spi=open(spifile,'w')
        spicmd='x70 = %s\n' %(numParts)
        spicmd+='do lb1 x10 = 1,x70\n' 
	spicmd+='FT\n'
	spicmd+='%s@{******x10}\n' %(stack)
	spicmd+='_1\n'
	spicmd+='MU\n' 
	spicmd+='_1\n'
	spicmd+='%s@{******x10}\n' %(ctfStack)
	spicmd+='_2\n'
	spicmd+='*\n'
	spicmd+='FT\n'
	spicmd+='_2\n'
	spicmd+='%s@{******x10}\n' %(outStack)
	spicmd+='lb1\n'
	runSpider(spicmd)

def runSpider(lines):
       spifile = "currentSpiderScript.spi"
       if os.path.isfile(spifile):
               os.remove(spifile)
       spi=open(spifile,'w')
       spi.write("MD\n")
       spi.write("TR OFF\n")
       spi.write("MD\n")
       spi.write("VB OFF\n")
       spi.write("MD\n")
       spi.write("SET MP\n")
       spi.write("(0)\n")
       spi.write("\n")
       spi.write(lines)

       spi.write("\nEN D\n")
       spi.close()
       spicmd = "spider spi @currentSpiderScript"
       spiout = subprocess.Popen(spicmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr.read()
       output = spiout.strip().split()
       if "ERROR" in output:
               print "Spider Error, check 'currentSpiderScript.spi'\n"
               sys.exit()
       #clean up
       os.remove(spifile)
       if os.path.isfile("LOG.spi"):
               os.remove("LOG.spi")
       resultf = glob.glob("results.spi.*")
       if resultf:
               for f in resultf:
                       os.remove(f)

def main(params):

	ctf = params['ctf']
	stack = params['stack']
	debug = params['debug']
	param= params['param']

	#Pixel size
	p13 = open(param,'r')
	pixel = 'pix'
	pixe = grep(pixel,p13)
	pi = pixe.split()
	pixSize = pi[2]	

	p2 = open(param,'r')
	#SNR
	n = 'snr'
	nr = grep(n,p2)
	rn = nr.split()
	amp = rn[2]

	#Box size
	p4 = open(param,'r')
	bxsz = 'boxSize'
	bxs = grep(bxsz,p4)
	bx = bxs.split()
	boxSize = bx[2]

	p5 = open(param,'r')
	#Number of particles
	nmpts = 'num_part'
	nmpt = grep(nmpts,p5)
	nmp = nmpt.split()
	numParts = nmp[2]

	#CS
	p6 = open(param,'r')
	cs1 = 'cs'
	cs2 = grep(cs1,p6)
	cs3 = cs2.split()
	Cs = cs3[2]
	

	#Convert stack into SPIDER format

	if debug is True:
		print 'proc2d %s %s.spi spider' %(stack,stack[:-4])	

	cmd = 'proc2d %s %s.spi spider' %(stack,stack[:-4])
	subprocess.Popen(cmd,shell=True).wait()

	#Convert parameter file into spider format

	f1 = open(ctf,'r')
	out = '%s.spi' %(ctf[:-4])
	if debug is True:
		print out
	o1 = open(out,'w')
	i = 1
	for line in f1:

		l = line.split()
		if debug is True:

			print line
			print l

		df1 = l[0]
		df2 = l[1]
		astig = l[2]

		o1.write('%s\t3\t%s\t%s\t%s\n' %(str(i),df1,df2,astig))
		i = i + 1

	o1.close()
	f1.close()

	#Run spider script to phaseflip particles
	#boxSize = str(92)
	#numParts = str(103)
	#Cs = str(6.2)
	#pixSize = (4.36)
	#amp = str(0.20)

	ctfStack = '%s_phase.spi' %(stack[:-4])

	flip(ctfStack[:-4],ctf[:-4],boxSize,numParts,Cs,float(pixSize),amp)
	outStack = '%s_flipped.spi' %(stack[:-4])
	mult(ctfStack[:-4],stack[:-4],outStack[:-4],numParts)
	
	cmd = 'proc2d %s %s.img' %(outStack,outStack[:-4])
	subprocess.Popen(cmd,shell=True).wait()

	#clean up
	cmd = 'rm %s %s %s.spi %s_phase.spi' %(out,outStack,stack[:-4],stack[:-4])
	subprocess.Popen(cmd,shell=True).wait()

if __name__ == "__main__":     
        params=setupParserOptions()     
        main(params)
