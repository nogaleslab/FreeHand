#!/usr/bin/env python

from sys import *
import os
from optparse import OptionParser
import glob
import subprocess
from os import system
import sys
import optparse 
from EMAN2 import *

def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog -p <parameter filek>")
        parser.add_option("-p",dest="param",type="string",metavar="FILE",
                help="EMAN2 parameter file")

        options,args = parser.parse_args()

        if len(args) > 1:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 2:
                parser.print_help()
                sys.exit()

        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params

def Eman2Freali(az,alt,phi):

    t1 = Transform({"type":"eman","az":az,"alt":alt,"phi":phi,"mirror":False})

    #t_conv = Transform({"type":"eman","alt":31.717474411458415,"az":90,"phi":-90,"mirror":False})

    #t2 = t1*t_conv.inverse()

    d = t1.get_params("eman")

    psi = d["phi"]+90

    if psi >360:

        psi = psi-360

    theta= d["alt"]

    phi = d["az"]-90

    return psi,theta,phi

def main(params):
	parm=params['param']

	f=open(parm,'r')
	out = open("%s_freeHand"%(parm),'w')
	count=1
	count2=1
	count=1
	for line in f:
					
		l = line.split()
		
		parmPSI = float(l[0])
		parmTHETA = float(l[1])
		parmPHI = float(l[2])
		sx =(float(l[3]))
		sy =(float(l[4]))
		model = float(l[5])

		psi,theta,phi = Eman2Freali(parmPSI,parmTHETA,parmPHI)	
			
		out.write("%s 	%s	%s	%s	%s	%s\n"%(psi,theta,phi,sx,sy,model))

		
	f.close()
	out.close()

#Need this at the end for the parse commands
if __name__ == "__main__":
     params=setupParserOptions()
     main(params)
