#!/usr/bin/env python

import linecache
import sys

#Convert parameter file format with CTF info

untilt = sys.argv[1] 
ctf2 = sys.argv[2]
fout = '%s_format' %(ctf2[:-4])
o1 = open(fout,'a')

o1.write("C Frealign format parameter file created from Search_fspace parameter file\n")
o1.write("C\n")
o1.write("C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMax\n")

i = 1

tmp = open(ctf2,'r')
tot = len(tmp.readlines())

while i <= tot:
	
	t = i + 3

	param = linecache.getline(untilt,t)
	ctf = linecache.getline(ctf2,i)

	l1 = param.split()
	l2 = ctf.split()

	psi = float(l1[1])
	theta = float(l1[2])
	phi = float(l1[3])

	shiftx = float(l1[4])
	shifty = float(l1[5])

	mag = float(l1[6])
	film=float(l1[7])
	
	df1 = float(l2[0])
	df2 = float(l2[1])
	astig = float(l2[2])
	a = (l1[10])
	test = '%s' %(a[-1:])

	if test == '*':
		CC = 50
	else:
		CC = float(l1[11])

	o1.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.0f%6d%9.1f%9.1f%8.2f%7.2f\n" %(i,psi,theta,phi,shiftx,shifty,mag,film,df1,df2,astig,CC))

	i = i + 1

o1.write("C\n")
