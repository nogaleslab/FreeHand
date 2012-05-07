#!/usr/bin/env python

import linecache
import sys

#Convert parameter file format with CTF info

#make_freeHand_Param_spi.py [angular file].spi [shifts file].spi [select file].spi [ctf].par [magnification] [pixel size]

ang = sys.argv[1]
shift = sys.argv[2]
select = sys.argv[3]
c = sys.argv[4]
mag = float(sys.argv[5])
pix = float(sys.argv[6])

fout2 = '%s.txt' %(select[:-4])
fout = '%s_format.par' %(c[:-4])
f1 = open(select,'r')
o1 = open(fout,'a')
o2 = open(fout2,'w')

o1.write("C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMax\n")

count = 1

for line in f1:

	l = line.split()
	t = l[0]
	if t[:1] is ';':	
		
		continue
	
	bl = l[2]
	bl2 = bl.strip('.000')
	
	good = int(bl2)

	o2.write('%s\n' %((good-1)))
	
	good = good + 1		#B/c first list has comments usually
	
	ag = linecache.getline(ang,good)
	ag = ag.split()
	psi = float(ag[2])
	theta = float(ag[3])
	phi = float(ag[4])

	sh = linecache.getline(shift,good)
	sh = sh.split()
	shiftx = float(sh[2])
	shifty = float(sh[3])
	cc = float(sh[4])*100

	ctf2 = linecache.getline(c,(good-1))
	ctf = ctf2.split()
	df1 = float(ctf[0])
	df2 = float(ctf[1])
	astig = float(ctf[2])

	o1.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.0f%6d%9.1f%9.1f%8.2f%7.2f\n" %(count,psi,theta,phi,shiftx,shifty,mag,1,df1,df2,astig,cc))

	count = count + 1

o1.write("C\n")
