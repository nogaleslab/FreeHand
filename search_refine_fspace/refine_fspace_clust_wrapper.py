#!/usr/bin/env python

##################################
#	Refinement info		##
##################################

model = 'ctf_file_00_merge_40A.mrc'	#Model name (mrc format)
stack = 'stack8_00.mrc'	#stack name (mrc format)
par = 'ctf_file_00_merge.par' 	#parameter file

#Parameters
apix = 6.02			#pixel size (A/pix)
amp = 0.07			#Amplitude contrast
cs = 2.2			#Spherical aberration
kv = 120			#Accelerating voltage

#refine_fspace inputs
thetaR = 20			#Theta constraint during refinement
psiR = 10			#psi constraint 
shR = 12				#shift constraint 
ri = 169			#radius for searching (ang)
rmax1 = 250			#low resolution limit (angstroms)
rmax2 = 35			#high resolution limit (angstroms)

#Number of processors
tot = 29936			#number of particle to refine
incr = 749			# nprocs = tot/incr
queue = 'barcelona.q'               #queue name
name = 'michael'                #name of user

###################################
###################################

import subprocess
import sys
import time
import os


#Get current working directory
script = sys.argv[0]
cwd = '%s' %(script[:-31])

i = 1
while i <= tot:
	n = i + incr
	if n > tot:	
		n = tot	

	cmd = 'qsub -q %s %s/refine_fspace_v1_02_clust.csh %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' %(queue,cwd,str(apix),str(amp),str(cs),str(kv),model,str(ri),str(rmax1),str(rmax2),str(thetaR),str(psiR),str(shR),str(i),str(n),stack,par,cwd )
	print cmd 
	subprocess.Popen(cmd,shell=True)
	time.sleep(3)
	i = i + incr + 1


