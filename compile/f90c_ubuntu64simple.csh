#!/bin/tcsh -x
###################################################################
#
# compile fortran program as xxx.f90 into xxx.exe using image libraries
#
###################################################################
#echo " starting compile procedure"

set prog=$1
set oput=$2
if ($prog == "" ) then
   echo " No input file given"
   echo " Must use : compile input_file"
   goto exit
endif
#
# if the input file (or with default extention .for)  exists?
#
set n=`find . -name "$prog" -print`
if ($#n == "0") then
   set prog1=$prog:t
   set prog2=$prog1:r
   set prog3=$prog2.f90
   set m=`find . -name "$prog3" -print`
   if ($#m == "0") then
      echo " xp: '$prog3' does not exsit"
      goto exit
   else
      set prog=$prog3
   endif
endif
#
# usage 2: naming output file with extention .exe
#
if ($oput == "") then
   set prog1=$prog:t
   set prog2=$prog1:r
   set oput=$prog2.exe
else
   set oput1=$oput:t
   set oput2=$oput1:r
   set oput=$oput2.com
endif
#
# The call for GNU's gfortran compiler
(time gfortran -O -w -o $oput $prog \
   /opt/qb3/mrc-2010/lib/imlib2010.a   \
   /opt/qb3/mrc-2010/lib/misclib.a   \
   /opt/qb3/mrc-2010/lib/genlib.a    \
    -lfftw3 -lm)
#
exit:
