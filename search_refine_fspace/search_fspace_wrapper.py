#!/usr/bin/env python

##################################
#	Refinement info		##
##################################

model = 'vols_2b_mr_002.mrc'	#Model name (mrc format)
stack = 'start.mrc'		#stack name (mrc format)
par = 'listCTFvalues.par_format2' #parameter file

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
tot = 100			#number of particle to refine
incr = 13			# nprocs = tot/incr

###################################
###################################

import subprocess

i = 1
while i <= tot:
	n = i + incr
	if n > tot:	
		n = maxi	

	cmd = 'search_fspace_v1_02.csh %f %f %f %f %f %f %s %f %f %f %f %f %f %s %s' %(i,n,apix,amp,cs,kv,model,ncheck,psi,sh,ri,rmax1,rmax2,stack,par)
	print cmd
	subprocess.Popen(cmd,shell=True)
	i = i + incr + 1

