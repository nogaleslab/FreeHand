#!/usr/bin/env python

import optparse
from sys import *
import os,sys,re
from optparse import OptionParser
import glob
import subprocess
from os import system
import linecache
import numpy as np
import matplotlib.pyplot as plt
import Image
import pylab 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#=========================
def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog -u <untilted stack> -m <model> -p <parameter file>")
        parser.add_option("-u",dest="untilted",type="string",metavar="FILE",
                help="untilted stack (white particles in IMAGIC format)")
	parser.add_option("-t",dest="tilted",type="string",metavar="FILE",
                help="tilted stack (black particles in IMAGIC format)")
        parser.add_option("-m",dest="model",type="string",metavar="FILE",
                help="3D model(s) for alignment (Single SPIDER volume or multi-volume HDF file)")
        parser.add_option("-p",dest="param",type="string", metavar="FILE",
                help="Parameter file with refinement info (free_param.par)")
        parser.add_option("-c",dest="ctf",type="string",metavar="FILE",
		help="CTF-information file for tilted particles")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 2:
                parser.print_help()
		sys.exit()
        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params

#=========================
def checkConflicts(params):
	if not params['untilted']:
		print "\nWarning: no untilted stack specified\n"
	elif not os.path.exists(params['untilted']):
		print "\nError: stack file '%s' does not exist\n" % params['untilted']
		sys.exit()
        if not params['tilted']:
                print "\nWarning: no tilted stack specified\n"
        elif not os.path.exists(params['tilted']):
                print "\nError: stack file '%s' does not exist\n" % params['tilted']
                sys.exit()        
	if not params['model']:
                print "\nWarning: no model specified\n"
        elif not os.path.exists(params['model']):
                print "\nError: model file '%s' does not exist\n" % params['model']
                sys.exit()
	if not params['param']:
		print "\nError: no free_param.par file specified"
		sys.exit()
	if not os.path.isfile(params['param']):
		print "\nError: free_param.par file does not exist\n" 
		sys.exit()
	if not os.path.isfile(params['ctf']):
                print "\nNo CTF-information specified for tilted stack; using 2 um as default\n"
                sys.exit()
	if os.path.exists('refine_eman2'):
		print "\nOutput folder already exists, remove refine_eman2 and start again\n"
		sys.exit()
	s = params['untilted']
	prep = '%s_prep.img' %(s[:-4])
	if os.path.exists(prep):
		os.remove(prep)
		os.remove('%s.hed' %(prep[:-4]))
	if os.path.exists('start.hdf'):
		os.remove('start.hdf')

#========================
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#========================
def getEMANPath():        
        ### get the imagicroot directory        
        emanpath = subprocess.Popen("env | grep EMAN2DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()        

        if emanpath:                
                emanpath = emanpath.replace("EMAN2DIR=","")                
        if os.path.exists(emanpath):                        
                return emanpath        
        print "EMAN2 was not found, make sure eman2/2.05 is in your path"        
        sys.exit()

#===========================
def getMYAMI():
	myamipath = subprocess.Popen("env | grep MYAMI", shell=True, stdout=subprocess.PIPE).stdout.read().strip()

        if myamipath:
                myamipath = myamipath.replace("MYAMI=","")
        if os.path.exists(myamipath):
                return myamipath
        print "myami was not found, make sure myami is loaded"
        sys.exit()

#========================
def getOPENMPIPath():
        ### get the openmpi directory        
        openpath = subprocess.Popen("env | grep MPIHOME", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        test = openpath.find('imagic')
	if test >= 0:
                print "OPENMPI is not loaded, make sure it is in your path"
                sys.exit()

	if test is None:

	        if openpath:
        	        openpath = openpath.replace("MPIHOME=","")
	        if os.path.exists(openpath):
        	        return openpath
	        print "OPENMPI is not loaded, make sure it is in your path"
        	sys.exit()

#==========================
def getCCP4Path():        
        ### get the openmpi directory        
        ccp4path = subprocess.Popen("env | grep CCP4_PATH", shell=True, stdout=subprocess.PIPE).stdout.read().strip()        

        if ccp4path:                
                ccp4path = ccp4path.replace("CCP4_PATH=","")      
        if os.path.exists(ccp4path):                        
                return ccp4path        
        print "ccp4 is not loaded, make sure it is in your path"        
        sys.exit()

#========================
def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            return match.string

#========================
def Eman2Freali(az,alt,phi):

    t1 = Transform({"type":"eman","az":az,"alt":alt,"phi":phi,"mirror":False})
    d = t1.get_params("eman")   
    psi = d["phi"]+90   
    if psi >360:        
        psi = psi-360
    theta= d["alt"]
    phi = d["az"]-90
    return psi,theta,phi

#========================
def align(params,cwd):
	debug = params['debug']
	param = params['param']
	untilt = params['untilted']
	model = params['model']

	#Get parameter info: angular step
	p = open(param,'r')
	a = 'angular' 
	angl = grep(a,p)
	aL = angl.split()
	ang = aL[2]
	
	#Shift
	s = 'shift'
	p = open(param,'r')
	shi = grep(s,p)
	sh = shi.split()
	sx = sh[2]

	#Pixel size
	p13 = open(param,'r')
	pixel = 'pix'
	pixe = grep(pixel,p13)
	pi = pixe.split()
	pix = pi[2]

	#Radius
	r = 'radius'
	p = open(param,'r')
	radiu = grep(r,p)
	radi = radiu.split()
	rad = radi[2]

	p2 = open(param,'r')
	#SNR
	n = 'snr'
	nr = grep(n,p2)
	rn = nr.split()
	snr = rn[2]

	p3 = open(param,'r')
	#ts
	ts1 = 'ts'
	ts2 = grep(ts1,p3)
	ts3 = ts2.split()
	ts = ts3[2]

	#Box size
	p4 = open(param,'r')
	bxsz = 'boxSize'
	bxs = grep(bxsz,p4)
	bx = bxs.split()
	box = bx[2]

	p5 = open(param,'r')
	#Number of particles
	nmpts = 'num_part'
	nmpt = grep(nmpts,p5)
	nmp = nmpt.split()
	tot = nmp[2]

	#CS
	p6 = open(param,'r')
	cs1 = 'cs'
	cs2 = grep(cs1,p6)
	cs3 = cs2.split()
	cs = cs3[2]

	#Accelerating voltage
	p7 = open(param,'r')
	v1 = 'volt'
	v2 = grep(v1,p7)
	v3 = v2.split()
	volt = v3[2]

	#Free hand angular search
	p8 = open(param,'r')
	fs1 = 'freeHand_ang_search'
	fs2 = grep(fs1,p8)
	fs3 = fs2.split()
	angSearch = fs3[2]
	
	#Free hand Low resolution  
	p9 = open(param,'r')
	mr1 = 'min_res'
	mr2 = grep(mr1,p9)
	mr3 = mr2.split()
	min_res = mr3[2]

	#Free hand Max resolution  
	p10 = open(param,'r')
	mr4 = 'min_res'
	mr5 = grep(mr4,p10)
	mr6 = mr5.split()
	max_res = mr6[2]

	#Free hand first particle
	p11 = open(param,'r')
	fir1 = 'first'
	fir2 = grep(fir1,p11)
	fir3 = fir2.split()
	first = fir3[2]

	#Free hand last particle
	p12 = open(param,'r')
	ls1 = 'last'
	ls2 = grep(ls1,p12)
	ls3 = ls2.split()
	last = ls3[2]

	#Free hand Max resolution  
	p10 = open(param,'r')
	mr4 = 'max_res'
	mr5 = grep(mr4,p10)
	mr6 = mr5.split()
	max_res = mr6[2]

        #Free hand increment  
        p13 = open(param,'r')
        inc1 = 'incr'
        inc2 = grep(inc1,p13)
        inc3 = inc2.split()
        incr = inc3[2]

        #Free hand increment  
        p14 = open(param,'r')
        m1 = 'mag'
        m2 = grep(m1,p14)
        m3 = m2.split()
        mag = m3[2]

        #Free hand increment  
        p15 = open(param,'r')
        m1 = 'num_mod'
        m2 = grep(m1,p15)
        m3 = m2.split()
        num_mod = int(m3[2])

	p17 = open(param,'r')
        pp1 = 'cutoff'
        pp2 = grep(pp1,p17)
        pp3 = pp2.split()
        cutoff = pp3[2]

        p18 = open(param,'r')
        pp1 = 'calc'
        pp2 = grep(pp1,p17)
        pp3 = pp2.split()
        calc = pp3[2]

		
	#Prepare stack for EMAN2 refinement 
        print '\n'
        print 'Converting stack into EMAN2 format'
        print '\n'

	#Filter particles to specified resolution limits
                
	cmd = 'proc2d %s %s_prep.img apix=%s hp=%s lp=%s' %(untilt,untilt[:-4],pix,min_res,max_res)
        subprocess.Popen(cmd,shell=True).wait()

       	cmd = '%s/up_head.py %s_prep.img %s' %(cwd,untilt[:-4],pix)
        subprocess.Popen(cmd,shell=True).wait()

        #Run refinement
        print '\n'                
	print 'Running EMAN2 refinement'                
	print '         Angular step = %s' %(ang)                
	print '         Shift range = %s' %(sx)                
	print '         Shift step size (ts)  = %s' %(ts)                
	print '         Pixel Size = %s' %(pix)                
	print '         Radius = %s' %(rad)                
	print '         SNR = %s' %(snr)
	print '	        CC_cut = %s' %(cutoff)                
	print '\n'                

	if num_mod == 1:

		cmd = 'mpirun -np 8 %s/refine.py start.hdf %s refine_eman2 --ou=%s --rs=1 --xr=%s --ts=%s --delta=%s --snr=%s --center=0 --maxit=1 --ref_a=S --sym=c1 --cutoff=%s --MPI --full_output' %(cwd,model,rad,sx,ts,ang,snr,cutoff)
		
		if debug is True:
			print cmd
               	subprocess.Popen(cmd,shell=True).wait()
	else:
	
		cmd = 'mpirun -np 8 %s/refine.py start.hdf %s refine_eman2 --ou=%s --rs=1 --xr=%s --ts=%s --delta=%s --snr=%s --center=0 --maxit=1 --ref_a=S --sym=c1 --cutoff=%s --MPI --full_output --sort' %(cwd,model,rad,sx,ts,ang,snr,cutoff)

                if debug is True:
                        print cmd
		subprocess.Popen(cmd,shell=True).wait()		
	
	#Clean up:
	cmd = 'rm logfile* start.hdf %s_prep.*' %(untilt[:-4])
	subprocess.Popen(cmd,shell=True).wait()

def eman2_sort(paramout,tilt,ctf,num_mod,debug):

	if debug is True:
		print 'eman2_sort():'
		print '		paramout = %s; tilt=%s; ctf=%s; num_mod=%s; debug=%s' %(paramout,tilt,ctf,num_mod,debug)

	#Sort particles by model(s)	
	if int(num_mod) == 1:
	
		if debug is True:
			print 'num_mod == 1'
       		param=open(paramout,'r')
       		count=1
       		text='%s_%02d.txt' %(tilt[:-4],0)
		c_o = '%s_model00.par' %(ctf[:-4])
		o1 = open(c_o,'w')
		y_o = '%s_model00' %(paramout)
		y1 = open(y_o,'w')

       		text=open(text,'w')

       		for line in param:
                	l=line.split()
                 	member=float(l[5])
			if debug is True:
				print l
                 	if member == 999:

                        	text.write("%s\n" %(count-1))
				
				c = linecache.getline(ctf,count)				
				y1.write('%s %s' %(str(count),line))				
				o1.write('%s' %(c))

 		  	count=count+1

      		text.close()
       		param.close()
	        cmd="proc2d %s %s_%02d.img list=%s_%02d.txt" %(tilt,tilt[:-4],0,tilt[:-4],0)
       		subprocess.Popen(cmd,shell=True).wait()

	else:

		for n in range(0,int(num_mod)):
			param=open(paramout,'r')
			c_o = '%s_model%02d.par' %(ctf[:-4],n)
                	o1 = open(c_o,'w')
			count=1
                	y_o = '%s_model%02d' %(paramout,n)
                	y1 = open(y_o,'w')
			text='%s_%02d.txt' %(tilt[:-4],n)	

			text=open(text,'w')

			for line in param:

				l=line.split()
				member=float(l[5])
	
				if member == n:

					text.write("%s\n" %(count-1))
					c = linecache.getline(ctf,count)
					y1.write('%s' %(line))
                         	        o1.write('%s' %(c))

				count=count+1
			text.close()
			param.close()
			cmd="proc2d %s %s_%02d.img list=%s_%02d.txt " %(tilt,tilt[:-4],n,tilt[:-4],n)
			subprocess.Popen(cmd,shell=True).wait()

def eman2_angConv(paramout,num_mod,ctf,mag,model,tilt,debug):
	mod_count = 0
	
	while mod_count < int(num_mod):

		print 'Working on model %s' %(mod_count)
			
		print '\n'                
		print 'Converting files into free-hand format'                
		print '\n'                

		parm='%s_model%02d' %(paramout,mod_count)
		if debug is True:
			print 'parm = %s' %(parm)

		f=open(parm,'r')
		out = open("%s_freeHand"%(parm),'w')
		count=1
		count2=1
		count=1
		for line in f:
					
			l = line.split()
			if debug is True:
				print l		
			parmPSI = float(l[1])
			parmTHETA = float(l[2])
			parmPHI = float(l[3])
			sx =(float(l[4]))
			sy =(float(l[5]))
			model1 = float(l[6])
			#Convert euler angles from EMAN2 to FREALIGN/SPIDER convention
			if debug is True:
				print 'parmPSI = %s	parmTHETA = %s	parmPHI = %s	' %(parmPSI,parmTHETA,parmPHI)
				
			psi,theta,phi = Eman2Freali(parmPSI,parmTHETA,parmPHI)	
			out.write("%s 	%s	%s	%s	%s	%s\n"%(psi,theta,phi,sx,sy,model1))
		
		f.close()
		out.close()

		makeFH_eman2('%s_freeHand' %(parm),'%s_model%02d.par' %(ctf[:-4],int(mod_count)),mag,1,debug)

		eman2_mods(num_mod,model,mod_count,debug)

		im_to_mrc('%s_%02d.img' %(tilt[:-4],mod_count),debug)

		mod_count = mod_count + 1

#=================
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

#============
def eman2_mods(num_mod,model,mod_count,debug):

        #Convert model from HDF to MRC                
	if debug is True:
		print num_mod
		print model
		print mod_count
	
        if int(num_mod) > 1:

                cmd = 'e2proc3d.py --first=%s --last=%s %s %s_%03d.mrc' %(model,model[:-4],mod_count)
        	if debug is True:
			print cmd
	        subprocess.Popen(cmd,shell=True).wait()

	else:

                cmd = 'proc3d %s %s_%03d.mrc' %(model,model[:-4],int(mod_count))
		if debug is True:
			print cmd
                subprocess.Popen(cmd,shell=True).wait()

#==================
def makeFH_eman2(f,c,mag,div,debug):

        #Convert parameter file format with CTF info
        f1 = open(f,'r')
        fout = '%s_format.par' %(f[:-4])
        o1 = open(fout,'a')
        if debug is True:
                print 'c = %s' %(c)
        o1.write("C Frealign format parameter file created from Search_fspace parameter file\n")
        o1.write("C\n")
        o1.write("C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMax\n")

        count = 1

        for line in f1:

                l = line.split()

                if debug is True:
                        print line
                psi = float(l[0])
                theta = float(l[1])
                phi = float(l[2])

                shiftx = float(l[3])/float(div)
                shifty = float(l[4])/float(div)

                ctf2 = linecache.getline(c,count)
                ctf = ctf2.split()
                df1 = float(ctf[0])
                df2 = float(ctf[1])
                astig = float(ctf[2])

                o1.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.0f%6d%9.1f%9.1f%8.2f%7.2f\n" %(count,psi,theta,phi,shiftx,shifty,float(mag),1,df1,df2,astig,50))

                count = count + 1

        o1.write("C\n")

#=========
def eman2_conv(params,cwd):
	
	param = params['param']

	#Get parameter info: number of models
	p = open(param,'r')
	a = 'num_mod' 
	angl = grep(a,p)
	aL = angl.split()
	num_mod = aL[2]
	
	#Get parameter info: mag
        p = open(param,'r')
        a = 'mag'
        angl = grep(a,p)
        aL = angl.split()
        mag = aL[2]

	p = open(param,'r')
        a = 'ang'
        angl = grep(a,p)
        aL = angl.split()
        ang = aL[2]

	tilt = params['tilted']
	ctf = params['ctf']
	debug = params['debug']
	model = params['model']
	paramout = 'refine_eman2/paramout_%03d_00' %(float(ang))

	#Sort particles based up membership to model(s)
	eman2_sort(paramout,tilt,ctf,num_mod,debug)
	
	#Convert euler angles, model, and particles from EMAN2 to FREALIGN for each model
	eman2_angConv(paramout,num_mod,ctf,mag,model,tilt,debug)

	#Clean up
	mod = 0
	while mod < int(num_mod):
		
		cmd = 'rm %s_model%02d %s_model%02d_freeHand %s_model%02d.par %s_%02d.img %s_%02d.hed %s_%02d.txt' %(paramout,mod,paramout,mod,ctf[:-4],mod,tilt[:-4],mod,tilt[:-4],mod,tilt[:-4],mod)
		if debug is True:
			print cmd
		subprocess.Popen(cmd,shell=True).wait()
		
		mod = mod + 1

def fastFree_run(params,cwd,mod):
	debug = params['debug']
	param = params['param']
	stack = params['tilted']
	model = params['model']

	#Pixel size
	p13 = open(param,'r')
	pixel = 'pix'
	pixe = grep(pixel,p13)
	pi = pixe.split()
	pix = pi[2]

	#Radius
	r = 'radius'
	p = open(param,'r')
	radiu = grep(r,p)
	radi = radiu.split()
	rad = radi[2]

	p2 = open(param,'r')
	#SNR
	n = 'snr'
	nr = grep(n,p2)
	rn = nr.split()
	snr = rn[2]

	#CS
	p6 = open(param,'r')
	cs1 = 'cs'
	cs2 = grep(cs1,p6)
	cs3 = cs2.split()
	cs = cs3[2]

	#Accelerating voltage
	p7 = open(param,'r')
	v1 = 'volt'
	v2 = grep(v1,p7)
	v3 = v2.split()
	volt = v3[2]

	#Free hand angular search
	p8 = open(param,'r')
	fs1 = 'freeHand_ang_search'
	fs2 = grep(fs1,p8)
	fs3 = fs2.split()
	angSearch = fs3[2]
	
	#Free hand Low resolution  
	p9 = open(param,'r')
	mr1 = 'min_res'
	mr2 = grep(mr1,p9)
	mr3 = mr2.split()
	min_res = mr3[2]

	#Free hand Max resolution  
	p10 = open(param,'r')
	mr4 = 'min_res'
	mr5 = grep(mr4,p10)
	mr6 = mr5.split()
	max_res = mr6[2]

	#Free hand first particle
	p11 = open(param,'r')
	fir1 = 'first'
	fir2 = grep(fir1,p11)
	fir3 = fir2.split()
	first = fir3[2]

	#Free hand last particle
	p12 = open(param,'r')
	ls1 = 'last'
	ls2 = grep(ls1,p12)
	ls3 = ls2.split()
	last = ls3[2]

	#Free hand Max resolution  
	p10 = open(param,'r')
	mr4 = 'max_res'
	mr5 = grep(mr4,p10)
	mr6 = mr5.split()
	max_res = mr6[2]

        #Free hand increment  
        p13 = open(param,'r')
        inc1 = 'incr'
        inc2 = grep(inc1,p13)
        inc3 = inc2.split()
        procs = inc3[2]

        #Free hand increment  
        p14 = open(param,'r')
        m1 = 'mag'
        m2 = grep(m1,p14)
        m3 = m2.split()
        mag = m3[2]

        #Free hand increment  
        p15 = open(param,'r')
        m1 = 'num_mod'
        m2 = grep(m1,p15)
        m3 = m2.split()
        num_mod = int(m3[2])

	p17 = open(param,'r')
        pp1 = 'cutoff'
        pp2 = grep(pp1,p17)
        pp3 = pp2.split()
        cutoff = pp3[2]

        p18 = open(param,'r')
        pp1 = 'calc'
        pp2 = grep(pp1,p17)
        pp3 = pp2.split()
        calc = pp3[2]

	p18 = open(param,'r')
        pp1 = 'angular'
        pp2 = grep(pp1,p18)
        pp3 = pp2.split()
        ang = pp3[2]

        #Number of particles

	ctf = 'refine_eman2/paramout_%03d_00_model%02d_free_format.par' %(float(ang),mod)

        y = open(ctf,'r')
        tot = len(y.readlines())
        tot = tot - 4

        
	#Run Free-Hand test                
	
	#Create ctf param info

	numParts = 0 
	dfCount = 0
	df = 0
	counter = 1
	ctfInfo = open(ctf)
	for line in ctfInfo:
		l = line.split()
	        if l[0] == 'C':
        		continue
		if debug is True:
			print l
        	df1 = float(l[8])
	        df2 = float(l[9])
        	astig = float(l[10])

		if counter > 1:
			if counter < tot:
			        if df != df1:
					try:
						 ctfEXE
					except NameError:
						ctfEXE = None
					if ctfEXE is None:
						ctfEXE = '%s,%s,1,%s,%s,%s,1\n' %(dfCount,mag,df1,df2,astig)
					else:
						ctfEXE += '%s,%s,1,%s,%s,%s,1\n' %(dfCount,mag,df1,df2,astig)

					df = df1
					dfCount = 0
			if counter == tot:
				dfCount = dfCount + 1
                	        ctfEXE += '%s,%s,1,%s,%s,%s,0\n' %(dfCount,mag,df1,df2,astig)

		if counter == 1:
			df = df1 
	
		dfCount = dfCount + 1
		counter = counter + 1

	i = 1
	        
	iteration = 1 

	incr = str(round(tot/int(procs))+1)

	#Setup inputs for free-hand executable
	cmd = 'cp %s/fastfreehand_v1_01.exe .' %(cwd)					
	subprocess.Popen(cmd,shell=True).wait()
	exe = 'fastfreehand_v1_01.exe\n'
	p1 = '%s,%s,%s,%s\n' %(pix,snr,cs,volt)					
	p2 = '1,0\n'
	p3 = '%s_%02d.mrc\n' %(stack[:-4],mod)
	p4 = '%s_%03d.mrc\n' %(model[:-4],mod)
	p5 = '%s\n' %(ctf)
	#p6 = plots_.mrc (see below)
	p7 = '-%s,%s,-%s,%s\n' %(angSearch,angSearch,angSearch,angSearch)
	p8 = '%s,%s,%s\n' %(min_res,max_res,str(float(pix)*float(rad)))
	#p9 = first, last (see below)
	p10 = '%s\n' %(calc)

        while i < int(tot): 
               	last = str(i + float(incr)-1)  
                last = last[:-2] 
                if i == 1:
        	        first = str(i) 
                else:
                        first = str(i)
	                first = first[:-2]
                if float(last) > int(tot): 
                        incr = int(incr[:-2]) - (int(last)- int(tot))
                        last = str(tot)
	
		p6 = 'model%02d_plots_CC_v101_%02d.mrc\n' %(mod,iteration)
		p9 = '%s,%s\n' %(first,last)

		ff_cmd ='#!/bin/csh\n'
		ff_cmd +='fastfreehand_v1_01.exe << eot\n'
		ff_cmd +=p1
		ff_cmd +=p2
		ff_cmd +=p3
		ff_cmd +=p4
                ff_cmd +=p5
                ff_cmd +=p6
                ff_cmd +=p7
                ff_cmd +=p8
                ff_cmd +=p9
                ff_cmd +=p10
                ff_cmd +=ctfEXE
		ff_cmd +='eot\n'
		ff_cmd +='touch iteration%01d_finished\n' %(iteration)

		tmp = open('tmp%01d.csh'%(iteration),'w')
		tmp.write(ff_cmd)
		tmp.close()

		cmd = 'chmod +x tmp%01d.csh' %(iteration)
		subprocess.Popen(cmd,shell=True).wait()
	
		cmd = './tmp%01d.csh' %(iteration)
		subprocess.Popen(cmd,shell=True)

		i = i + float(incr)
               	iteration = iteration + 1
	iteration = iteration - 1
	return iteration

def fastFree(params,cwd):

	param = params['param']

	#Free hand increment  
        p15 = open(param,'r')
        m1 = 'num_mod'
        m2 = grep(m1,p15)
        m3 = m2.split()
        num_mod = int(m3[2])

	mod = 1
	while mod <= num_mod:
		mod = mod -1 
		iteration = fastFree_run(params,cwd,mod)
		wait(params,iteration)
		mod = mod + 2
#===========	
def wait(params,iteration):

        param = params['param']
        debug = params['debug']

        i = 1

        while i<= iteration:

                test = os.path.isfile('iteration%01d_finished'%(i))

                if test is False:
                        time.sleep(5)
			
                if test is True:
                        i = i + 1

        if debug is True:
                print 'Free-hand test completed for all particles'

        #Clean up:
        cmd = 'rm iteration?_finished tmp* fastfreehand_v1_01.exe'
        subprocess.Popen(cmd,shell=True)

#==============
def findPeak(stack,peakfile):
	if os.path.exists(peakfile):
        	os.remove(peakfile)
        out = open(peakfile,'w')
        stackRead = mrc.read(stack)
        number, ydim, xdim = stackRead.shape
        for i in range(number):
        	image = stackRead[i,:,:]
                output = peakfinder.findPixelPeak(image, guess=None, limit=None, lpf=None)
                coord = output['pixel peak']
                out.write('%d   %s      %s\n' %(i,coord[0],coord[1]))
	return out

#==============
def averageStack(stackfile,avgfile):
               
                a = mrc.read(stackfile)
                a = np.sum(a,axis=0)
                a = (a-a.min())/(a.max()-a.min())
                mrc.write(a,avgfile)
                return avgfile

#==============
def convertEPS_to_PNG(image,out):

	im = Image.open(image)
        im.rotate(-90).save(out)

#===============
def scatter(data,lim,tilt,include):

        tiltX = tilt[0]
        tiltY = tilt[1]

        loadData = np.loadtxt(data)
        x = loadData[:,2]
        y = loadData[:,1]

        print x
	print y
	#Center peak vales at (0,0)
        centx = np.subtract(x,lim)
        centy = np.subtract(y,lim)

        #Calculate distance of peaks from expected angle
        dist = []
        for i in xrange(len(loadData[:,1])):
	        rx = centx[i]
        	ry = centy[i]
		print rx
                distance = math.sqrt(((rx - tiltX)*(rx - tiltX)) + ((ry - tiltY)*(ry - tiltY))/2)
                dist.append(distance)

        newDist = sorted(dist)

        numReadLines = round((float(include)/100)*len(loadData[:,1]))

        includeRadius = []
        for j in xrange(numReadLines):
        	includeRadius = newDist[j]
        #Create function for plotting circle
        theta = np.linspace(0,2*math.pi)

        incRadx = includeRadius*pylab.cos(theta)
        incRady = includeRadius*pylab.sin(theta)

        incRadx = pylab.add(incRadx,tiltX)
        incRady = pylab.add(incRady,tiltY)

        #Set radii for concentric circles       
        rad1 = 10
        rad2 = 20
        rad3 = 30
        rad4 = 40
        rad5 = 50
        rad6 = 60

        #Create x,y coords for concentric cricles
        cx1 = rad1*pylab.cos(theta)
        cy1 = rad1*pylab.sin(theta)
        cx2 = rad2*pylab.cos(theta)
        cy2 = rad2*pylab.sin(theta)
        cx3 = rad3*pylab.cos(theta)
        cy3 = rad3*pylab.sin(theta)
        cx4 = rad4*pylab.cos(theta)
        cy4 = rad4*pylab.sin(theta)
        cx5 = rad5*pylab.cos(theta)
        cy5 = rad5*pylab.sin(theta)
        cx6 = rad6*pylab.cos(theta)
        cy6 = rad6*pylab.sin(theta)

        #Create axes
        line1 = np.linspace(-60,60,100)

        #Create zeros array
        line2 = []
        i = 1
        while i <= 100:
                line2.append('0')
		i = i + 1
        fig = plt.figure(1)

        scatter = plt.subplot(111,aspect='equal')
        majorLocator = MultipleLocator(10)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(5)
        scatter.set_xlabel('Tilt direction (degrees)',fontsize=15)
        scatter.set_ylabel('Tilt direction (degrees)',fontsize=15)
        scatter.set_title('%d'%(include) + '% ' + 'of particles are within %d degrees of expected angle'%(round(includeRadius)))
        scatter.plot(cx1,cy1,c = 'k')
        scatter.plot(cx2,cy2,c = 'k')
        scatter.plot(cx3,cy3,c = 'k')
        scatter.plot(cx4,cy4,c = 'k')
        scatter.plot(cx5,cy5,c = 'k')
        scatter.plot(cx6,cy6,c = 'k')
        scatter.plot(cx5,cy5,c = 'k')
        scatter.plot(incRadx,incRady,c = 'r')
        scatter.plot(line2,line1,c='k')
        scatter.plot(line1,line2,c='k')
        scatter.scatter(centx,centy,marker='+',c='k',edgecolor='k',s=55)
        majorLocator = MultipleLocator(10)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(5)
        scatter.xaxis.set_major_locator(majorLocator)
        scatter.xaxis.set_major_formatter(majorFormatter)
        scatter.xaxis.set_minor_locator(minorLocator)
        scatter.yaxis.set_major_locator(majorLocator)
        scatter.yaxis.set_major_formatter(majorFormatter)
        scatter.yaxis.set_minor_locator(minorLocator)
        plt.xlim(-lim,lim)
        plt.ylim(-lim,lim)
	outFILE = '%s.png' %(data[:-4])
        plt.savefig(outFILE,dpi=150,format='png')
#===========
def plotFH(params,ccp4_path,cwd):

	param = params['param']
	debug = params['debug']
	model = params['model']
	stack = params['tilted']

	#Free hand increment  
        p15 = open(param,'r')
        m1 = 'num_mod'
        m2 = grep(m1,p15)
        m3 = m2.split()
        num_mod = int(m3[2])

        #Free hand angular search
        p8 = open(param,'r')
        fs1 = 'freeHand_ang_search'
        fs2 = grep(fs1,p8)
        fs3 = fs2.split()
        angSearch = fs3[2]

        p18 = open(param,'r')
        fs1 = 'calc'
        fs2 = grep(fs1,p18)
        fs3 = fs2.split()
        calc = fs3[2]

	p18 = open(param,'r')
        fs1 = 'tilt_ang'
        fs2 = grep(fs1,p18)
        fs3 = fs2.split()
        tilt_ang = int(fs3[2])

	tiltCenter = [tilt_ang,0]
	includedPercentTilt = 40

	mod = 1

	while mod <= num_mod:
		mod = mod - 1
		#Merge stacks:
		m = 0

		num = len(glob.glob('model%02d_plots*.mrc'%(mod)))

		i = 1
		while i <= int(num):

			cmd = 'e2proc2d.py model%02d_plots_CC_v101_%02d.mrc model%02d_plots_CC_v101_%02d.img --threed2twod' %(mod,i,mod,i)
			if debug is True:
				print cmd
			subprocess.Popen(cmd,shell=True).wait()

			cmd = 'rm model%02d_plots_CC_v101_%02d.mrc' %(mod,i)
			subprocess.Popen(cmd,shell=True).wait()
	
			cmd = 'proc2d model%02d_plots_CC_v101_%02d.img model%02d_plots_CC_v101_merge.img' %(mod,i,mod)
			if debug is True:
				print cmd
			subprocess.Popen(cmd,shell=True).wait()

			i = i + 1

		cmd = 'e2proc2d.py model%02d_plots_CC_v101_merge.img model%02d_plots_CC_v101_merge.mrc --twod2threed' %(mod,mod)
		if debug is True:
			print cmd
		subprocess.Popen(cmd,shell=True).wait()

		cmd = 'cp %s/totsumstack.exe .' %(cwd)	
		subprocess.Popen(cmd,shell=True).wait()

		totsum = '#!/bin/csh\n'
		totsum += 'totsumstack.exe << eot\n'
		totsum += 'model%02d_plots_CC_v101_merge.mrc\n' %(mod) 
		totsum += 'model%02d_averageplot_CC_v101.mrc\n' %(mod)
		totsum += 'eot\n'

		tmp = open('tmp.csh','w')
		tmp.write(totsum)
		tmp.close()

		if debug is True:
		        print totsum

		cmd = 'chmod +x tmp.csh' 
		if debug is True:
		        print cmd
		subprocess.Popen(cmd,shell=True)

		cmd = './tmp.csh' 
		subprocess.Popen(cmd,shell=True).wait()

		cmd = 'rm totsumstack.exe tmp.csh '
		subprocess.Popen(cmd,shell=True).wait()
	
#		avgfile = averageStack('model%02d_plots_CC_v101_merge.mrc' %(mod),'model%02d_averageplot_CC_v101.mrc'%(mod)) 

		if calc is 'C':

			line1 = (float(angSearch)*2)/5
			line = line1/2
			
			npo = '#!/bin/csh\n'
			npo += 'rm -f z.plot\n'
			npo += 'rm -f plot84.ps\n'
			npo += '%s/bin/npo mapin model%02d_averageplot_CC_v101.mrc plot z.plot << eof\n' %(ccp4_path,mod)
			npo += 'NOTITLE\n'
			npo += 'MAP SCALE 1 INVERT\n'
			npo += '# For CCC\n'
			npo += 'CONTRS 0.0 to 1 by 0.002\n'
			npo += 'LIMITS 0 %s 0 %s 0 0\n' %(str(float(angSearch)*2),str(float(angSearch)*2))
			npo += 'SECTNS 0 0 1\n'
			npo += 'GRID  5 5\n'
			npo += 'GRID U DASHED 1.0 0.2 0 EVERY %s FULL\n' %(line)
			npo += 'GRID V DASHED 1.0 0.2 0 EVERY %s FULL\n' %(line)
			npo += 'PLOT Y\n'
			npo += 'eof\n'
			npo += '\n'
			npo += '%s/bin/pltdev -log -dev ps -abs -pen c -xp 3.1 -yp 3.1 -lan -i z.plot -o model%02d_average_frehand_CC.ps\n' %(ccp4_path,mod)

			tmp = open('tmp.csh','w')
			tmp.write(npo)
			tmp.close()

			cmd = 'chmod +x tmp.csh'
			subprocess.Popen(cmd,shell=True)

			cmd = './tmp.csh'
			subprocess.Popen(cmd,shell=True).wait()

			cmd = 'rm tmp.csh z.plot'
			#subprocess.Popen(cmd,shell=True).wait()

			findPeak('model%02d_plots_CC_v101_merge.mrc' %(mod),'model%02d_peaks.txt' %(mod))
			#scatter('model%02d_peaks.txt' %(mod),angSearch,tiltCenter,includedPercentTilt)
			convertEPS_to_PNG('model%02d_average_frehand_CC.ps'%(mod),'model%02d_average_frehand_CC.png' %(mod))
			os.remove('model%02d_average_frehand_CC.ps' %(mod))

			cmd = 'rm -r model%02d_plots_CC_v101_merge.* model%02d_plots_CC_v101_??.* %s_%03d.mrc %s_%02d.mrc' %(mod,mod,model[:-4],mod,stack[:-4],mod)
			subprocess.Popen(cmd,shell=True).wait()
	
		if calc is 'P':

		        line1 = (float(angSearch)*2)/5
		        line = line1/2
		        
		        npo = '#!/bin/csh\n'
		        npo += 'rm -f z.plot\n'
		        npo += 'rm -f plot84.ps\n'
		        npo += '%s/bin/npo mapin model%02d_averageplot_CC_v101.mrc plot z.plot << eof\n' %(ccp4_path,mod)
		        npo += 'NOTITLE\n'
		        npo += 'MAP SCALE 1 INVERT\n'
		        npo += '# For Pres\n'
		        npo += 'CONTRS 77. to 86. by .3\n'
		        npo += 'LIMITS 0 %s 0 %s 0 0\n' %(str(float(angSearch)*2),str(float(angSearch)*2))
		        npo += 'SECTNS 0 0 1\n'
		        npo += 'GRID  5 5\n'
		        npo += 'GRID U DASHED 1.0 0.2 0 EVERY %s FULL\n' %(line)
		        npo += 'GRID V DASHED 1.0 0.2 0 EVERY %s FULL\n' %(line)
		        npo += 'PLOT Y\n'
		        npo += 'eof\n'
			npo += '\n'
			npo += '%s/bin/pltdev -log -dev ps -abs -pen c -xp 3.1 -yp 3.1 -lan -i z.plot -o model%02d_average_frehand_CC.ps\nn' %(ccp4_path,mod)

		        tmp = open('tmp.csh','w')
		        tmp.write(npo)
		        tmp.close()

		        cmd = 'chmod +x tmp.csh'
		        subprocess.Popen(cmd,shell=True)

		        cmd = './tmp.csh'
		        subprocess.Popen(cmd,shell=True).wait()

		        cmd = 'rm tmp.csh z.plot'
		        subprocess.Popen(cmd,shell=True).wait()

			cmd = 'proc2d model%02d_plots_CC_v101_merge.img model%02d_plots_CC_v101_merge.img invert inplace' %(mod,mod)
			subprocess.Popen(cmd,shell=True).wait()

		        cmd = 'e2proc2d.py model%02d_plots_CC_v101_merge.img model%02d_plots_CC_v101_merge.spi' %(mod,mod)
		        subprocess.Popen(cmd,shell=True).wait()

		        tot = EMUtil.get_image_count('model00_plots_CC_v101_merge.img')  
		        n = int(angSearch)+1
		        stack1 = 'model%02d_plots_CC_v101_merge.spi'%(mod)
		        peak(stack1,tot,n)

			cmd = 'rm -r model%02d_plots_CC_v101_merge.* model%02d_plots_CC_v101_??.* %s_%03d.mrc %s_%02d.mrc' %(mod,mod,model[:-4],mod,stack[:-4],mod)
		        subprocess.Popen(cmd,shell=True).wait()

		mod = mod + 2

if __name__ == "__main__":     
	getMYAMI()
	from pyami import peakfinder,mrc
	getEMANPath()             
	getOPENMPIPath()
	ccp4 = getCCP4Path()
	from EMAN2 import *     
	from sparx  import *     
	params=setupParserOptions()     
	checkConflicts(params)
	cwd = '/archive/michael/lib/freeHand'
	align(params,cwd)
	eman2_conv(params,cwd)
	fastFree(params,cwd)
	plotFH(params,ccp4,cwd)
	cmd = 'rm -r refine_eman2'
	subprocess.Popen(cmd,shell=True).wait()
