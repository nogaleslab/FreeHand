#!/usr/bin/env python

#Output filenames
paramOUT1 = 'parameter_00.par'	#param filename for untilted particles
paramOUT2 = 'parameter_01.par'	#param filename for tilted particles
stack1 = 'stack00'		#stack name for untilted particles
stack2 = 'stack01'		#stack name for tilted particles

#Inputs
shrink = 4			#Binning factor for final particle stack
new = 64			#Box size for binned particles
scale = 4			#Binning factor used for micrograph from which the particles were picked

#CTFTILT inputs
parm3 = "2.2,120.0,0.2,99993,15,2\n" # !CS[mm],HT[kV],AmpCnst,XMAG,DStep[um]
parm4 = "128,400.0,8.0,2000.0,30000.0,500.0,30,5\n" #!Box,ResMin[A],ResMax[A],dFMin[A],dFMax[A],FStep

#############		Script		###############

import subprocess
import sys
import re
import os
import glob

p1 = open(paramOUT1,'wa')
p2 = open(paramOUT2,'wa')
cmd = "/opt/qb3/frealign-9.02/bin/ctftilt.exe"

#syntax: grep(regexp_string,list_of_strings_to_search)
def grep(string,list):
        expr = re.compile(string)
        return filter(expr.search,list)

list = glob.glob('*en_00.mrc')

for file in list: 
	tmp1 = re.sub("_00_","_01_",file)
       	tiltname = re.sub("en_00","en_01",tmp1) 
       	#Check if both files exist
	if os.path.isfile(file) & os.path.isfile(tiltname):
               	a = subprocess.Popen(cmd, -1, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		[o,e] = a.communicate('%s\n'%(file) + '%s.diag\n'%(file) + parm3 + parm4)
		out = grep("Final Values", o.split("\n"))
		out2 = out[0]
		out3 = out2.split()
		
		fname = file.strip('.mrc')
       		b = open('%s.box'%(fname))
        	tot = len(b.readlines())
        	i = 1
        
        	while i <= tot:
	
			p1.write('%s		%s		%s		%s\n' %(out3[0],out3[1],out3[2],out3[4]))
			i = i + 1

		cmd2 = 'batchboxer input=%s.mrc dbbox=%s.box scale=%s output=%s.img' %(fname,fname,scale,stack1)
		subprocess.Popen(cmd2,shell=True).wait()

		a = subprocess.Popen(cmd, -1, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                [o,e] = a.communicate('%s\n'%(tiltname) + '%s.diag\n'%(tiltname) + parm3 + parm4)
                out = grep("Final Values", o.split("\n"))
                out2 = out[0]
                out3 = out2.split()

                fname = tiltname.strip('.mrc')
                b = open('%s.box'%(fname))
                tot = len(b.readlines())
                i = 1

                while i <= tot:

                        p2.write('%s    	%s      	%s      	%s\n' %(out3[0],out3[1],out3[2],out3[4]))
                        i = i + 1

                cmd2 = 'batchboxer input=%s.mrc dbbox=%s.box scale=%s output=%s.img' %(fname,fname,scale,stack2)
                subprocess.Popen(cmd2,shell=True).wait()

cmd3 = 'proc2d %s.img %s_dc4.img meanshrink=%s' %(stack1,stack1,shrink)
subprocess.Popen(cmd3,shell=True).wait()

cmd4 = 'proc2d %s_dc4.img %s_dc4_%03d.img clip=%s edgenorm=0,1' %(stack1,stack1,float(new),new)
subprocess.Popen(cmd4,shell=True).wait()

cmd3 = 'proc2d %s.img %s_dc4.img meanshrink=%s' %(stack2,stack2,shrink)
subprocess.Popen(cmd3,shell=True).wait()

cmd4 = 'proc2d %s_dc4.img %s_dc4_%03d.img clip=%s edgenorm=0,1' %(stack2,stack2,float(new),new)
subprocess.Popen(cmd4,shell=True).wait()
