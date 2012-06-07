#!/usr/bin/env python

# Author: Richard J. Hall 04/2011 rjhall@berkeley.edu 

import os
from global_def import *
from optparse import OptionParser
import sys

from EMAN2  import *

def main():
	arglist = []
        for arg in sys.argv:
		arglist.append( arg )
	        progname = os.path.basename(arglist[0])
	        usage = progname + " raw_images apix "
	        parser = OptionParser(usage,version="1.0")
	        (options, args) = parser.parse_args(arglist[1:])
        if len(args) != 2:
                print "usage: " + usage
                print "Please run '" + progname + " -h' for detailed options"
	else:
		a = EMData()
		fname = args[0]
		apix = float(args[1])
		imn = EMUtil.get_image_count(fname) 
		for i in xrange(imn):
			a.read_image(fname,i)
			a.set_attr_dict({'active':1})
			t2 = Transform({"type":"spider","phi":0,"theta":0,"psi":0})
		        a.set_attr("xform.projection", t2)
			a.set_attr("apix_x",apix )
			a.write_image("start.hdf",i)

if __name__ == "__main__":
        main()

