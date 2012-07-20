#!/usr/bin/env python

import linecache
import sys

#Convert parameter file format with CTF info

f = sys.argv[1] 
f1 = open(f,'r')
fout = '%s_format' %(f[:-4])
o1 = open(fout,'a')
mag = float(sys.argv[2])

o1.write("C Frealign format parameter file created from Search_fspace parameter file\n")
o1.write("C\n")
o1.write("C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMax\n")

count = 1

for line in f1:

	l = line.split()
	
	psi = float(l[1])
	theta = float(l[2])
	phi = float(l[3])

	shiftx = float(l[4])
	shifty = float(l[5])

	film=float(l[7])
	
	df1 = float(l[8])
	df2 = float(l[9])
	astig = float(l[10])

	o1.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.0f%6d%9.1f%9.1f%8.2f%7.2f\n" %(count,psi,theta,phi,shiftx,shifty,mag,film,df1,df2,astig,0))

	count = count + 1

o1.write("C\n")
