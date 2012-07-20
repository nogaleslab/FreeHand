#!/usr/bin/env python

import subprocess
import sys
import time

l = sys.argv[1]
tot = float(sys.argv[2])
incr = float(sys.argv[3])

new = '%smerge.par' %(l)

o1 = open(new,'a')

o1.write("C Frealign format parameter file created from Search_fspace parameter file\n")
o1.write("C\n")
o1.write("C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMax\n")

i = 1

while i <= tot:

	n = i + incr
	if n > tot:
		n = tot

	if i == 1:
		p1 = '1'
	if i > 1:
		p1 = str(i)
		p1 = p1[:-2]
	p2 = str(n)

	f1 = '%s%s_%s' %(l,p1,p2[:-2])

	q1=open(f1,'r')

	for qline in q1:
		t = qline.split()
		t = t[0]	
		if t[:1] is 'C':
			continue
		o1.write('%s' %(qline))

	q1.close()

	i = i + incr + 1

q1.close()

