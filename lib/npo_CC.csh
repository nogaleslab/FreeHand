\rm -f z.plot
\rm -f plot84.ps
/public/ccp4/bin/npo \
MAPIN plotsum_CC.mrc \
PLOT z.plot \
<< 'END-npo'
NOTITLE
MAP SCALE 1 INVERT
# For CCC
CONTRS 0.0 to 1 by 0.002
# For Pres
#CONTRS 77. to 86. by .3
#MODE BELOW 60.0 DASHED 1 0.10 0
LIMITS 0 90 0 90 0 0
SECTNS 0 0 1
#GRID 5 5 
#GRID 45 45
GRID  5 5
GRID U DASHED 1.0 0.2 0 EVERY 9 FULL
GRID V DASHED 1.0 0.2 0 EVERY 9 FULL
PLOT Y
END-npo

/public/ccp4/bin/pltdev -log -dev ps -abs -pen c -xp 3.1 -yp 3.1 -i z.plot -o frehand_CC.ps
