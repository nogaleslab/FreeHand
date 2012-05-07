#!/bin/csh -f

set stack = $1
set model = $2
set sx = $3
set ang = $4
set radius = $5
set snr = $6
set ts = $7
set cutoff = $8

mpirun -np 8 /opt/qb3/eman-2.0-rc3/contrib/nogales/refine.py $stack $model refine_eman2 --ou=$radius --rs=1 --xr=$sx --ts=$ts --delta=$ang --snr=$snr --center='0' --maxit=1 --ref_a=S --sym=c1 --cutoff=$cutoff --MPI --debug --full_output

echo mpirun -np 8 /opt/qb3/eman-2.0-rc3/contrib/nogales/refine.py $stack $model refine_eman2 --ou=$radius --rs=1 --xr=$sx --ts=$ts --delta=$ang --snr=$snr --center='0' --maxit=1 --ref_a=S --sym=c1 --cutoff=$cutoff --MPI --debug --full_output
