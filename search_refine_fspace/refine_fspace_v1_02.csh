#!/bin/csh -x 

#READ (5,*)  psize,wgh,cs,akv
#READ (5,*)  ctfexppart,ctfexpmod
#READ (5,*)  inmap      ! input map
#READ (5,*)  ri
#READ (5,*)  rmax1,rmax2
#READ (5,*)  thetaphibound,psibound,shiftbound ! deg,deg,pix
#READ (5,*)  ifirst,ilast,nstacks
#READ (5,*)  instack
#READ (5,*)  inparfile
#READ (5,*)  outparfile

set apix = $1
set amp = $2
set cs = $3
set volt = $4
set model = $5
set ri = $6
set rmax1 = $7
set rmax2 = $8
set theta = $9
set psi = $10
set sh = $11
set first = $12
set last = $13
set stack = $14
set par = $15
set p = $16

time $p/refine_fspace_v1_02.exe << eot
$apix,$amp,$cs,$volt
1,2
$model
$ri
$rmax1,$rmax2
$theta,$psi,$sh
$first,$last,1
$stack
$par
${par}_refined_${first}_${last}.par 
eot
