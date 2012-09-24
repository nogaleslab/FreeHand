#!/bin/csh 

# Name of job
#$ -N refine_fspace

# Set parallel environment; set number of processors
#$ -pe smp 1

# Max walltime for this job (2 hrs)
##$ -l h_rt=02:00:00

# Merge the standard out and standard error to one file
##$ -j y

# Run job through csh shell
#$ -S /bin/csh

# use current working directory
#$ -cwd

# The following is for reporting only. It is not really needed
# to run the job. It will show up in your output file.
#

echo "Job starting `date`"
echo "Current working directory: $cwd"
echo "Got $NSLOTS processors."

# The job



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
${par:r}_refined_${first}_${last}
eot
