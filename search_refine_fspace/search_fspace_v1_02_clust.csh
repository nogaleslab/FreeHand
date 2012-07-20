#!/bin/csh 

# Name of job
#$ -N search_fspace

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


# psize (e.g 2.4),wgh (0.07),cs (2.0), akv (200.0) 
# ctfexppart,ctfexpmod
# inmap      ! input map
# ncheck (1000),psi_step (5) ,shiftmax (pixels: 5),ri (170)
# rmax1 (300),rmax2 (20)
# ifirst (1) ,ilast (19867), nstacks (1)
# instack(1)  ! input particle stack
# inputparfile(1)
# outputparfile(1) 

set first = $1
set last = $2
set apix = $3
set amp = $4
set cs = $5
set kv = $6
set model = $7
set ncheck = $8
set psi = $9
set sh = $10
set ri = $11
set rmax1 = $12
set rmax2 = $13
set stack = $14
set par = $15


time search_fspace_v1_02.exe << eot
$apix,$amp,$cs,$kv
1,2
$model
$ncheck,$psi,$sh,$ri
$rmax1,$rmax2
$1,$2,1
$stack
$par
${par}_${1}_${2}
eot

