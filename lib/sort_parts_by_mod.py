#!/usr/bin/env python 

import sys
import subprocess
from optparse import OptionParser
import optparse
import linecache

def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog -f <stack> -p <parameter file> --num=<FLOAT>")
        parser.add_option("-f",dest="stack",type="string",metavar="FILE",
                help="IMAGIC particle stack")
        parser.add_option("-p",dest="param",type="string",metavar="FILE",
                help="EMAN2 output parameter file")
	parser.add_option("-c",dest="ctf",type="string",metavar="FILE",
                help="Per particle CTF-information")
        parser.add_option("--num",dest="num",type="int", metavar="INT",
                help="number of models used in refinement")
        parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
	options,args = parser.parse_args()
	

        if len(args) > 4:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 4:
                parser.print_help()
                sys.exit()
	params={}        
	for i in parser.option_list: 
        	if isinstance(i.dest,str):
	                params[i.dest] = getattr(options,i.dest)        
	return params
def main(params):

	stack=params['stack'] 
	numMods=int(params['num'])
	new= stack[:-4]
	ctf = params['ctf']

	if numMods == 1:

		debug = params['debug']
       		param=open(params['param'],'r')
       		count=1
       		text='%s_%02d.txt' %(new,0)
		c_o = '%s_model00.par' %(ctf[:-4])
		o1 = open(c_o,'w')
		y_o = '%s_model00' %(params['param'])
		y1 = open(y_o,'w')

       		text=open(text,'w')

       		for line in param:
                	l=line.split()
                 	member=float(l[5])
			if debug is True:
				print l
                 	if member == 999:

                        	text.write("%s\n" %(count-1))
				
				c = linecache.getline(ctf,count)				
				y1.write('%s' %(line))				
				o1.write('%s' %(c))

 		  	count=count+1

      		text.close()
       		param.close()
	        cmd="proc2d %s %s_%02d.img list=%s_%02d.txt" %(stack,new,0,new,0)
       		subprocess.Popen(cmd,shell=True).wait()
	else:
		for n in range(0,numMods):
			param=open(params['param'],'r')
			c_o = '%s_model%02d.par' %(ctf[:-4],n)
                	o1 = open(c_o,'w')
			count=1
                	y_o = '%s_model%02d' %(params['param'],n)
                	y1 = open(y_o,'w')
			text='%s_%02d.txt' %(new,n)	

			text=open(text,'w')

			for line in param:

				l=line.split()
				member=float(l[5])
	
				if member == n:

					text.write("%s\n" %(count-1))
					c = linecache.getline(ctf,count)
					y1.write('%s' %(line))
                         	        o1.write('%s' %(c))

				count=count+1
			text.close()
			param.close()
			cmd="proc2d %s %s_%02d.img list=%s_%02d.txt " %(stack,new,n,new,n)
			subprocess.Popen(cmd,shell=True).wait()

if __name__ == "__main__":
     params=setupParserOptions()
     main(params)
