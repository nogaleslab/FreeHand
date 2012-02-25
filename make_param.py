#!/usr/bin/env python

import sys

o = open(sys.argv[1],'w')
num = float(sys.argv[2])

i = 1

df1 = 10639
df2 = 10255
asti = -59.310
mag = 49000
zero = 0
while i <= num: 

	o.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.0f%6d%9.1f%9.1f%8.2f%7.2f%7.2f\n" %(i,zero,zero,zero,zero,zero,mag,1,df1,df2,asti,0,0))

	i = i + 1
