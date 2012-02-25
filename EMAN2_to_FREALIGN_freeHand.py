#!/usr/bin/env python

import optparse
from sys import *
import os,sys,re
from optparse import OptionParser
import glob
import subprocess
from os import system
import linecache

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -f <stack> -p <parameter> -c <ctf> --num=<float> --mag=<float>")
	parser.add_option("-f",dest="stack",type="string",metavar="FILE",
		help="raw, IMAGIC particle stack (black particles)")
	parser.add_option("-p",dest="param",type="string",metavar="FILE",
		help="EMAN2 output parameter file")
	parser.add_option("-c",dest="ctf",type="string",metavar="FILE",
		help="per-particle CTF information file from ctftilt.py")
	parser.add_option("--num",dest="num",type="int", metavar="INT",
		help="number of models used in refinement")
	parser.add_option("--mag",dest="mag",type="float", metavar="FLOAT",
		help="actual magnification of images")
        parser.add_option("--box",dest="box",type="float", metavar="FLOAT",
                help="Box size")
	parser.add_option("-s", action="store_true",dest="spider",default=False,
		help="Flag if your original stack before converting to HDF format was a SPIDER stack OR if you used XMIPP normalization during APPION processing")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 5:
		parser.print_help()
		
	params={}
	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params

def getIMAGICPath():
	### get the imagicroot directory
	impath = subprocess.Popen("env | grep IMAGIC_ROOT", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        imagicpath = impath.replace("IMAGIC_ROOT=","")
        if imagicpath != '/opt/qb3/imagic-070813':
			print "imagic/070813 was not found, make sure it is in your path"
        		sys.exit()

def getEMANPath():        
	### get the imagicroot directory        
	emanpath = subprocess.Popen("env | grep EMAN2DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()        
	
	if emanpath:                
		emanpath = emanpath.replace("EMAN2DIR=","")                
	if os.path.exists(emanpath):                        
		return emanpath        
	print "EMAN2 was not found, make sure it is in your path"        
	sys.exit()

def getImagicVersion(imagicroot):
        ### get IMAGIC version from the "version_######S" file in
        ### the imagicroot directory, return as an int
        versionstr=glob.glob(os.path.join(imagicroot,"version_*"))
        if versionstr:
                v = re.search('\d\d\d\d\d\d',versionstr[0]).group(0)
                return int(v)
        else:
                print "Could not get version number from imagic root directory"
		sys.exit()

def Eman2Freali(az,alt,phi):

    t1 = Transform({"type":"eman","az":az,"alt":alt,"phi":phi,"mirror":False})

    #t_conv = Transform({"type":"eman","alt":31.717474411458415,"az":90,"phi":-90,"mirror":False})

    #t2 = t1*t_conv.inverse()

    d = t1.get_params("eman")

    psi = d["phi"]+90

    if psi >360:

        psi = psi-360

    theta= d["alt"]

    phi = d["az"]-90

    return psi,theta,phi

def main(params):
			
		parm=params['param']
		numMods = params['num']
		mag = params['mag']
		stack = params['stack']
		debug = params['debug']
		box = str(params['box'])

		if numMods == 1:
			f=open(parm,'r')
			out = open("%s_01_frealign"%(parm),'w')
		  	count=1
		        text2='%s_%02d.txt' %(parm,1)
      			text=open(text2,'w')
			count2=1
			
			print "\n"
			print "Calculating euler angle conversion..."
			print "\n"
				
			for line in f:
					
				l = line.split()
	 	
				parmPSI = float(l[0])
				parmTHETA = float(l[1])
				parmPHI = float(l[2])
				sx =(float(l[3]))
				sy =(float(l[4]))
				model = float(l[5])

				psi,theta,phi = Eman2Freali(parmPSI,parmTHETA,parmPHI)	

				if model == 999:
					if debug is True:

						print 'Particle %s is included' %(count-1)					
	
					text.write("%s\n" %(count-1))

					ctf = linecache.getline(params['ctf'],count)
					if debug is True:

						print 'Reading line %s in ctf file' %(count)
						print ctf
					c = ctf.split()

					df1 = float(c[0])
					df2 = float(c[1])
					astig = float(c[2])					

					out.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.1f%6d%9.1f%9.1f%8.2f%7.2f%6.2f\n" %(count2,psi,theta,phi,sx,sy,mag,1,df1,df2,astig,0,0))
					count2 = count2 + 1

				count = count + 1

			f.close()
			out.close()
			text.close()

			new = stack[:-4]
		
			print "\n"
			print "Selecting particles..."
			print "\n"

			cmd="e2proc2d.py --list=%s %s %s_model1.img" %(text2,stack,new)
      			subprocess.Popen(cmd,shell=True).wait()

			print "\n"
			print "Normalizing selected particles..."
			print "\n"

			cmd="proc2d %s_model1.img %s_model1_norm.img norm=0,1" %(new,new)
			subprocess.Popen(cmd,shell=True).wait()

			if params['spider'] is True:

				print "\n"
				print "Flipping particle coordinates..."
				print "\n"

				cmd="proc2d %s_model1_norm.img %s_model1_norm.img inplace flip" %(new,new)
				subprocess.Popen(cmd,shell=True).wait()

			print "\n"
			print "Converting IMAGIC stack into FREALIGN particle stack..."
			print "\n"

			newO = open("%s_01_frealign"%(parm),'r')
			tot = len(newO.readlines())
			newO.close()
	
			cmd="~michael/BATCHLIB/freeHand/testim.b %s %s" %(box[:-2],tot)
			subprocess.Popen(cmd,shell=True).wait()

			cmd="rm %s_model1_norm.hed" %(new)
			subprocess.Popen(cmd,shell=True).wait()

			cmd="cp test.hed %s_model1_norm.hed" %(new)
			subprocess.Popen(cmd,shell=True).wait()

			cmd="~michael/BATCHLIB/freeHand/imagic_to_3Dmrc.b %s_model1_norm.img" %(new)
			subprocess.Popen(cmd,shell=True).wait()

		if numMods == 2:
			f=open(parm,'r')
			out0 = open("%s_00_frealign"%(parm),'w')
		  	out1= open("%s_01_frealign"%(parm),'w')
			count=1
		        text2='%s_%02d.txt' %(parm,0)
      			text3='%s_%02d.txt' %(parm,1)
			text=open(text2,'w')
			text1=open(text3,'w')			
			count2=1
			count3=1
			
			print "\n"
			print "Calculating euler angle conversion..."
			print "\n"
				
			for line in f:
					
				l = line.split()
	 	
				parmPSI = float(l[0])
				parmTHETA = float(l[1])
				parmPHI = float(l[2])
				sx =(float(l[3]))
				sy =(float(l[4]))
				model = float(l[5])

				psi,theta,phi = Eman2Freali(parmPSI,parmTHETA,parmPHI)	

				if model == 0:
					if debug is True:

						print 'Particle %s is included' %(count-1)					
	
					text.write("%s\n" %(count-1))

					ctf = linecache.getline(params['ctf'],count)
					if debug is True:

						print 'Reading line %s in ctf file' %(count)
						print ctf
					c = ctf.split()

					micro = float(c[7])
					df1 = float(c[8])
					df2 = float(c[9])
					astig = float(c[10])					

					out0.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.1f%6d%9.1f%9.1f%8.2f%7.2f%6.2f\n" %(count2,psi,theta,phi,sx,sy,mag,micro,df1,df2,astig,0,0))
					count2 = count2 + 1

				if model == 1:
					if debug is True:

						print 'Particle %s is included' %(count-1)					
	
					text1.write("%s\n" %(count-1))

					ctf = linecache.getline(params['ctf'],count)
					if debug is True:

						print 'Reading line %s in ctf file' %(count)
						print ctf
					c = ctf.split()

					micro = float(c[7])
					df1 = float(c[8])
					df2 = float(c[9])
					astig = float(c[10])					

					out1.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.1f%6d%9.1f%9.1f%8.2f%7.2f%6.2f\n" %(count3,psi,theta,phi,sx,sy,mag,micro,df1,df2,astig,0,0))
					count3 = count3 + 1

				count = count + 1

			f.close()
			out1.close()
			out0.close()
			text.close()
			text1.close()

			new = stack[:-4]
		
			print "\n"
			print "Selecting particles..."
			print "\n"

			cmd="e2proc2d.py --list=%s %s %s_model00.img" %(text2,stack,new)
      			subprocess.Popen(cmd,shell=True).wait()

			cmd="e2proc2d.py --list=%s %s %s_model01.img" %(text3,stack,new)
      			subprocess.Popen(cmd,shell=True).wait()

			print "\n"
			print "Normalizing selected particles..."
			print "\n"

			cmd="proc2d %s_model00.img %s_model00_norm.img norm=0,1" %(new,new)
			subprocess.Popen(cmd,shell=True).wait()
			
			cmd="proc2d %s_model01.img %s_model01_norm.img norm=0,1" %(new,new)
			subprocess.Popen(cmd,shell=True).wait()

			if params['spider'] is True:

				print "\n"
				print "Flipping particle coordinates..."
				print "\n"

				cmd="proc2d %s_model00_norm.img %s_model00_norm.img inplace flip" %(new,new)
				subprocess.Popen(cmd,shell=True).wait()

				cmd="proc2d %s_model01_norm.img %s_model01_norm.img inplace flip" %(new,new)
				subprocess.Popen(cmd,shell=True).wait()

			print "\n"
			print "Converting IMAGIC stack into FREALIGN particle stack..."
			print "\n"

			cmd="~michael/bin/imagic_to_frealign.b %s_model00_norm.img" %(new)
			subprocess.Popen(cmd,shell=True).wait()
			
			cmd="~michael/bin/imagic_to_frealign.b %s_model01_norm.img" %(new)
			subprocess.Popen(cmd,shell=True).wait()	
				
			

if __name__ == "__main__":
     imagicroot = getIMAGICPath()
     getEMANPath()
     from EMAN2 import *
     from sparx  import *
     params=setupParserOptions()
     main(params)
   
