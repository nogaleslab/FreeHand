#!/usr/bin/env python

import linecache
import sys

#Convert parameter file format with CTF info

f = sys.argv[1] 
f1 = open(f,'r')
fout = '%s_format' %(f)
o1 = open(fout,'a')

o1.write("C Frealign format parameter file created from Search_fspace parameter file\n")
o1.write("C\n")
o1.write("C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMax\n")

c = sys.argv[2] 
mag = float(sys.argv[3])

count = 1

for line in f1:

	l = line.split()
	
	psi = float(l[0])
	theta = float(l[1])
	phi = float(l[2])

	shiftx = float(l[3])
	shifty = float(l[4])

	ctf2 = linecache.getline(c,count)
	ctf = ctf2.split()
	df1 = float(ctf[0])
	df2 = float(ctf[1])
	astig = float(ctf[2])

	o1.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.0f%6d%9.1f%9.1f%8.2f%7.2f\n" %(count,psi,theta,phi,shiftx,shifty,mag,1,df1,df2,astig,50))

	count = count + 1

o1.write("C\n")
