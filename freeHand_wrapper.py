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
        parser.set_usage("%prog -u <untilted stack> -t <tilted stack> -m <model> -p <parameter file> -c <ctf info>")
        parser.add_option("-u",dest="untilted",type="string",metavar="FILE",
                help="untilted stack (white particles in IMAGIC format)")
        parser.add_option("-t",dest="tilted",type="string",metavar="FILE",
                help="tilted stack (black particles in IMAGIC format)")
        parser.add_option("-m",dest="model",type="string",metavar="FILE",
                help="3D model(s) for alignment (Single SPIDER volume or multi-volume HDF file)")
        parser.add_option("-p",dest="param",type="string", metavar="FILE",
                help="Parameter file with refinement info (free_param.par)")
        parser.add_option("-c",dest="ctf",type="string", metavar="FILE",
                help="Per-particle CTF info")
        parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 4:
                parser.print_help()
		sys.exit()
        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

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
        print "EMAN2 was not found, make sure eman2/2.0-rc3 is in your path"        
        sys.exit()

def getOPENMPIPath():        
        ### get the openmpi directory        
        openpath = subprocess.Popen("env | grep MPIHOME", shell=True, stdout=subprocess.PIPE).stdout.read().strip()        

        if openpath:                
                openpath = openpath.replace("MPIHOME=","")      
        if os.path.exists(openpath):                        
                return openpath        
        print "openmpi is not loaded, make sure it is in your path"        
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

def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            return match.string

def main(params):
	debug = params['debug']
	param = params['param']
	untilt = params['untilted']
	tilt = params['tilted']
	model = params['model']
	ctf = params['ctf']

	#Get current working directory
	script = sys.argv[0]
	cwd = '%s/lib' %(script[:-20])

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
                
	cmd = 'proc2d %s %s_prep.img apix=%s hp=%s lp=%s' %(untilt,untilt[:-4],pix,min_res,max_res)
        subprocess.Popen(cmd,shell=True).wait()

	if debug is True:
		print '%s/up_head.py %s_prep.img %s' %(cwd,untilt[:-4],pix)

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
		
		if debug is True:
			print '%s/run_freehand.sh start.hdf %s %s %s %s %s %s %s %s' %(cwd,model,str(sx),str(ang),str(rad),str(snr),str(ts),str(cutoff),cwd)

		cmd = '%s/run_freehand.sh start.hdf %s %s %s %s %s %s %s %s' %(cwd,model,str(sx),str(ang),str(rad),str(snr),str(ts),str(cutoff),cwd)
               	subprocess.Popen(cmd,shell=True).wait()
	else:
		if debug is True:
			print '%s/run_freehand_sort.sh start.hdf %s %s %s %s %s %s %s %s' %(cwd,model,sx,ang,rad,snr,ts,cutoff,cwd)
	
		cmd = '%s/run_freehand_sort.sh start.hdf %s %s %s %s %s %s %s %s' %(cwd,model,sx,ang,rad,snr,ts,cutoff,cwd)
		subprocess.Popen(cmd,shell=True).wait()		

	paramout = 'paramout_%03d_00' %(int(ang))
		
	if debug is True:
		print'%s/sort_parts_by_mod.py -f %s -p refine_eman2/%s -c %s --num=%s' %(cwd,tilt,paramout,ctf,num_mod)

	cmd = '%s/sort_parts_by_mod.py -f %s -p refine_eman2/%s -c %s --num=%s' %(cwd,tilt,paramout,ctf,num_mod)
	subprocess.Popen(cmd,shell=True).wait()

	mod_count = 0
	
	while mod_count < num_mod:

		print 'Working on model %s' %(mod_count)
			
		print '\n'                
		print 'Converting files into free-hand format'                
		print '\n'                
		paramout = 'paramout_%03d_00_model%02d' %(int(ang),mod_count)                
	
		if debug is True:
			print '%s/EMAN_to_FREALIGN_fromParam_freeHand.py -p refine_eman2/%s' %(cwd,paramout)

		cmd = '%s/EMAN_to_FREALIGN_fromParam_freeHand.py -p refine_eman2/%s' %(cwd,paramout)                
		subprocess.Popen(cmd,shell=True).wait()                
		
		#Convert parameter file format with CTF and angular info                
		cmd = '%s/make_freeHand_Param.py refine_eman2/%s_freeHand %s_model%02d.par %s' %(cwd,paramout,ctf[:-4],mod_count,mag)                
		subprocess.Popen(cmd,shell=True).wait()                

		#Convert model from HDF to MRC                
				
		if num_mod > 1:

			if debug is True:
				print '%s/e2proc3d_all.py %s %s' %(cwd,model,mod_count)

			cmd = '%s/e2proc3d_all.py %s %s' %(cwd,model,mod_count)
			subprocess.Popen(cmd,shell=True).wait()

			cmd = 'rm %s_%03d.spi' %(model[:-4],mod_count)
			subprocess.Popen(cmd,shell=True).wait()

			cmd = 'mv %s_%03d_mr.spi %s_%03d.spi' %(model[:-4],mod_count,model[:-4],mod_count)
			subprocess.Popen(cmd,shell=True).wait()

		else:
			if debug is True:
				print 'cp %s %s_%03d.spi' %(model,model[:-4],mod_count)

			cmd = 'cp %s %s_%03d.spi' %(model,model[:-4],mod_count)
			subprocess.Popen(cmd,shell=True).wait()

		if debug is True:
			print '%s/3D_spi_to_3D_mrc.b %s_%03d.spi' %(cwd,model[:-4],mod_count)

		cmd = '%s/3D_spi_to_3D_mrc.b %s_%03d.spi' %(cwd,model[:-4],mod_count)                
		subprocess.Popen(cmd,shell=True).wait()                

		tot = file_len('refine_eman2/%s' %(paramout)) 

		#Convert tilted particles to 3D-MRC format                

		if debug is True:
			print '%s/imagic_to_freeHand2.py -f %s_%02d.img --total=%s --box=%s' %(cwd,tilt[:-4],mod_count,tot,box)
			
		cmd = '%s/imagic_to_freeHand2.py -f %s_%02d.img --total=%s --box=%s' %(cwd,tilt[:-4],mod_count,tot,box)                
		subprocess.Popen(cmd,shell=True).wait()

               
		#Run Free-Hand test                
		info = linecache.getline('refine_eman2/%s_freeHand_format' %(paramout),4)                

		i = info.split()                

		#mag = i[6]                
		df1 = i[8]
	        df2 = i[9]
                astig = i[10]
    
	        i = 1
	        
		iteration = 1 

               	while i < int(tot): 
                	last = str(i + float(incr)-1)  
                        last = last[:-2] 
                        if i == 1:
                                 first = str(i) 
                        else:
                                 first = str(i)
	                         first = first[:-2]
                        if float(last) > int(tot): 
                                 incr = int(incr) - (int(last)- int(tot))
                                 last = str(tot)

			if debug is True:
				print '%s/fastFreeHand_wrapper.csh %s %s %s %s %s_%02d.mrc %s_%03d.mrc refine_eman2/%s_freeHand_format %s %s %s %s %s %s %s %s %s %s %s %s model%02d %s' %(cwd,pix,snr,cs,volt,tilt[:-4],mod_count,model[:-4],mod_count,paramout,angSearch,min_res,max_res,str(float(pix)*float(rad)),first,last,incr,mag,df1,df2,astig,str(iteration),mod_count,calc)
				

                        if debug is False:
				cmd= '%s/fastFreeHand_wrapper.csh %s %s %s %s %s_%02d.mrc %s_%03d.mrc refine_eman2/%s_freeHand_format %s %s %s %s %s %s %s %s %s %s %s %s model%02d %s' %(cwd,pix,snr,cs,volt,tilt[:-4],mod_count,model[:-4],mod_count,paramout,angSearch,min_res,max_res,str(float(pix)*float(rad)),first,last,incr,mag,df1,df2,astig,str(iteration),mod_count,calc)
                        	subprocess.Popen(cmd,shell=True)

       	                i = i + float(incr)
               	        iteration = iteration + 1


		mod_count = mod_count + 1

if __name__ == "__main__":     
	imagicroot = getIMAGICPath()     
	getEMANPath()             
	getOPENMPIPath()
	from EMAN2 import *     
	from sparx  import *     
	params=setupParserOptions()     
	main(params)


