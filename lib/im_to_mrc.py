#!/usr/bin/env 

from EMAN2 import *
from sparx import *
import sys

def im_to_mrc(stack,debug):

        #Convert tilted particles to 3D-MRC format                

        # get box size
        im=EMData.read_images(stack,[0])
        nx = im[0].get_xsize()
        del im
        nimg = EMUtil.get_image_count(stack)

        img = EMData(nx,nx,nimg)
        img.write_image(stack[:-4]+'.mrc')

        i = 0

        while i < nimg:
                d = EMData()
                d.read_image(stack, i)
                region = Region(0, 0, i, nx, nx, 1)
                d.write_image(stack[:-4]+".mrc",0,EMUtil.get_image_ext_type("mrc"), False, region, EMUtil.EMDataType.EM_FLOAT, True)
                i = i + 1

f = sys.argv[1]
debug = False

im_to_mrc(f,debug)
