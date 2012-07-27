#!/usr/bin/env python 

#This will convert all volumes within an HDF file into separate SPIDER volumes

#To run:

#./e2proc3d_all.py volume_file.hdf [volume number wanted]

import glob
import subprocess
import sys

f = sys.argv[1]
num=float(sys.argv[2])

new=f[:-4]

cmd="e2proc3d.py --first=%01d --last=%01d %s %s_%03d.hdf" %(num,num,f,new,num)
subprocess.Popen(cmd,shell=True).wait()
		
cmd="proc3d %s_%03d.hdf %s_%03d.img imagic" %(new,num,new,num)
subprocess.Popen(cmd,shell=True).wait()

cmd="/home/michael/BATCHLIB/freeHand/e2proc3d.b %s_%03d.img" %(new,num)
subprocess.Popen(cmd,shell=True).wait()

cmd="rm %s_%03d.hdf" %(new,num)
subprocess.Popen(cmd,shell=True).wait()

cmd="rm %s_%03d.hed %s_%03d.img" %(new,num,new,num)
subprocess.Popen(cmd,shell=True).wait()

cmd="~michael/BATCHLIB/freeHand/mirror_vols_spi.py %s_%03d.spi" %(new,num)
subprocess.Popen(cmd,shell=True).wait()

 
