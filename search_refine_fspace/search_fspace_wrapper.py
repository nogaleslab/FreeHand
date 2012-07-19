#!/usr/bin/env python

import subprocess

i = 1
maxi = 10000 
inr = 1250

while i <= maxi:

	n = i + inr

	if n > maxi:	
		n = maxi	

	cmd = 'search_fspace_v1_02.csh %s %s' %(i,n)
	print cmd
	subprocess.Popen(cmd,shell=True)
	i = i + inr + 1


