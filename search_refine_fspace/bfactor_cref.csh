#READ (5,*)  psize,molweight(kDa),scatfact,bfactor
#READ (5,*)  inmap      ! input map
#READ (5,*)  outmap
#READ (5,*)  infsc
time /home/jlr/programs/bfactor_cref/bfactor_cref.exe << eot
2.8,600,0.42,-1000
thirdmap.mrc
thirdmap_filt.mrc
1,fsc.txt
16
eot
