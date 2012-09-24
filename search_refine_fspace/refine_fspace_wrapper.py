#!/usr/bin/env python

##################################
#	Refinement info		##
##################################

model = 'iid_iia_iib_scp_4.mrc'	#Model name (mrc format)
stack = 'start00_01.mrc'	#stack name (mrc format)
par = 'parameter_00_sel_01_format_merge.par' 	#parameter file

#Parameters
apix = 6.02			#pixel size (A/pix)
amp = 0.07			#Amplitude contrast
cs = 2.2			#Spherical aberration
kv = 120			#Accelerating voltage

#refine_fspace inputs
thetaR = 60			#Theta constraint during refinement
psiR = 10			#psi constraint 
shR = 12				#shift constraint 
ri = 169			#radius for searching (ang)
rmax1 = 250			#low resolution limit (angstroms)
rmax2 = 30			#high resolution limit (angstroms)

#Number of processors
tot = 81			#number of particle to refine
incr = 11			# nprocs = tot/incr

###################################
###################################

import subprocess
import sys
import time
import os


#Get current working directory
script = sys.argv[0]
cwd = '%s' %(script[:-25])

i = 1
while i <= tot:
	n = i + incr
	if n > tot:	
		n = tot	

	cmd = '%s/refine_fspace_v1_02.csh %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' %(cwd,str(apix),str(amp),str(cs),str(kv),model,str(ri),str(rmax1),str(rmax2),str(thetaR),str(psiR),str(shR),str(i),str(n),stack,par,cwd )
	print cmd 
	subprocess.Popen(cmd,shell=True)
	i = i + incr + 1

i = i - 1 - incr

#Waiting loop:
time.sleep(10)
i = 1

while i <= tot:

	n = i + incr
	if n > tot:
		n = tot

	f = '%s_refined_%s_%s' %(par,i,n)
	fsize = 0
	while fsize == 0:
	
		test = os.path.exists(f)
		if test is False:
			print "Error, %s does not exist; must be bug in refine_fspace" %(f)
			sys.exit()
	
		fsize = os.path.getsize(f)
		time.sleep(1)

		if fsize > 0:
			i = i + incr + 1

tmp = str(float(incr))

cmd = '%s/combine_parfiles.py %s_refined_ %s %s' %(cwd,par,tot,tmp[:-2])
subprocess.Popen(cmd,shell=True).wait()	

cmd = 'rm %s_refined_*_*' %(par)
subprocess.Popen(cmd,shell=True).wait()

