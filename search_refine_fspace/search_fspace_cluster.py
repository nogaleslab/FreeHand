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
ncheck = 200			#ncheck, how many tries? don't know
psi = 5				#psi step for searching (?)
sh = 10				#shift search limit (pix)
ri = 35				#radius for searching (pix)
rmax1 = 300			#low resolution limit (angstroms)
rmax2 = 30			#high resolution limit (angstroms)

#Number of processors
tot = 30800			#number of particle to refine
incr = 963			# nprocs = tot/incr
queue = 'himem.q' 		#queue name
name = 'michael'		#name of user

###################################
###################################

import subprocess
import sys
import commands
import time

def qstat(name):
        q = commands.getoutput('qstat -u "*"')
        t1 = q.find(name)

        if t1 < 0:
                return 'False'

        t2 = t1 + len(name) + 12

        p = q[t1:t2]

        p1 = p.split()

        if p1[1] is 'qw':
                return 'qw'
        if p1[1] is 'r':
                return 'True'

i = 1
while i <= tot:
	n = i + incr
	if n > tot:	
		n = tot	
	
	cmd = 'qsub -q %s search_fspace_v1_02_clust.csh %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' %(queue,str(i),str(n),str(apix),str(amp),str(cs),str(kv),model,str(ncheck),str(psi),str(sh),str(ri),str(rmax1),str(rmax2),stack,par)
	print cmd
	subprocess.Popen(cmd,shell=True)
	
	time.sleep(8)
	i = i + incr + 1
