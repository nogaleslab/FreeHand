
#READ (5,*)  psize,wgh,cs,akv
#READ (5,*)  ctfexppart
#READ (5,*)  rmax2
#READ (5,*)  ifirst,ilast,nstacks
#READ (5,*)  outfile1,maskfile,threshargument (</>),thresh ! map
#DO i=1,nstacks
#   READ (5,*)  instack(i) ! particle stack
#   READ (5,*)  parfile(i) ! input Parameters
#ENDDO
time build_fspace_v1_01.exe << eot
6.02,0.07,2.2,120
1
30
1,10000,1
test_v1_01_3Dstack2.mrc,1,test3d.mrc,>,0.18
start_3Dstack2.mrc
listCTFvalues_format_merge.par
eot
# 10000
