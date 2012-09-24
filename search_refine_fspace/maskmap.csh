#READ (5,*) inmapfile
#READ (5,*) outmaskfile
#READ (5,*) outmapfile
#READ (5,*) psize (A/pix),rmax2 (angstroms),thresh (stddevs),length(pixels)
time /home/jlr/programs/maskmap/maskmap.exe << eot
3Dmap_V1Vo_3FLAG_apix_4.2_lp_21.mrc
mask.mrc
map.mrc
2.8,20,3,3
eot