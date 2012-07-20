#!/usr/bin/env python

import subprocess
import sys
import time

l = sys.argv[1]

cmd = 'ls %s* > %s' %(l,'outputtest.log')
subprocess.Popen(cmd,shell=True)

time.sleep(1)

f1 = open('outputtest.log','r')
new = '%smerge.par' %(l)

o1 = open(new,'a')

o1.write("C Frealign format parameter file created from Search_fspace parameter file\n")
o1.write("C\n")
o1.write("C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMax\n")

for line in f1:
	
	l = line.split()
	q = l[0]

	q1=open(q,'r')

	print 'Working on file %s' %(q)

	for qline in q1:
		t = qline.split()
		t = t[0]	
		if t[:1] is 'C':
			continue
		o1.write('%s' %(qline))

	q1.close()

f1.close()

cmd = 'rm outputtest.log'
subprocess.Popen(cmd,shell=True)
