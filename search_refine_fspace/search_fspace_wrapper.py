#!/usr/bin/env python

##################################
#	Refinement info		##
##################################

model = 'vols_2b_mr_002.mrc'	#Model name (mrc format)
stack = 'start.mrc'		#stack name (mrc format)
par = 'listCTFvalues_format' #parameter file

#Parameters
apix = 6.02			#pixel size (A/pix)
amp = 0.07			#Amplitude contrast
cs = 2.2			#Spherical aberration
kv = 120			#Accelerating voltage

#search_fspace inputs
ncheck = 10			#ncheck, how many tries? don't know
psi = 5				#psi step for searching (?)
sh = 10				#shift search limit (pix)
ri = 35				#radius for searching (pix)
rmax1 = 300			#low resolution limit (angstroms)
rmax2 = 30			#high resolution limit (angstroms)

#Number of processors
tot = 1000			#number of particle to refine
incr = 130			# nprocs = tot/incr

###################################
###################################

import subprocess
import sys

#Get current working directory
script = sys.argv[0]
cwd = '%s' %(script[:-25])

i = 1
while i <= tot:
	n = i + incr
	if n > tot:	
		n = tot	

	cmd = '%s/search_fspace_v1_02.csh %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' %(cwd,str(i),str(n),str(apix),str(amp),str(cs),str(kv),model,str(ncheck),str(psi),str(sh),str(ri),str(rmax1),str(rmax2),stack,par,cwd)
	print cmd 
	subprocess.Popen(cmd,shell=True)
	i = i + incr + 1

