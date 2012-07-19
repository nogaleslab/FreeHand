# psize (e.g 2.4),wgh (0.07),cs (2.0), akv (200.0) 
# ctfexppart,ctfexpmod
# inmap      ! input map
# ncheck (1000),psi_step (5) ,shiftmax (pixels: 5),ri (170)
# rmax1 (300),rmax2 (20)
# ifirst (1) ,ilast (19867), nstacks (1)
# instack(1)  ! input particle stack
# inputparfile(1)
# outputparfile(1) 
time search_fspace_v1_02.exe << eot
6.02,0.07,2.2,120
1,2
vols_2b_mr_002.mrc
100,5,10,35
300,30
$1,$2,1
start.mrc
listCTFvalues.par_format
listCTFvalues.par_format_out_${1}_${2}
eot
