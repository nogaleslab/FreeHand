#!/bin/csh

set limit = $1
set f = $2
set line = $3 

\rm -f z.plot
\rm -f plot84.ps
/opt/qb3/ccp4-6.1/bin/npo mapin ${f}_averageplot_CC_v101.mrc plot z.plot << eof
NOTITLE
MAP SCALE 1 INVERT
# For CCC
#CONTRS 0.0 to 1 by 0.002
# For Pres
CONTRS 77. to 86. by .3
#MODE BELOW 60.0 DASHED 1 0.10 0
LIMITS 0 $limit 0 $limit 0 0
SECTNS 0 0 1
#GRID 5 5 
#GRID 45 45
GRID  5 5
GRID U DASHED 1.0 0.2 0 EVERY $line FULL
GRID V DASHED 1.0 0.2 0 EVERY $line FULL
PLOT Y
eof

/opt/qb3/ccp4-6.1/bin/pltdev -log -dev ps -abs -pen c -xp 3.1 -yp 3.1 -lan -i z.plot -o ${f}_frehand_CC.ps
rm z.plot
