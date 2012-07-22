#!/usr/bin/env python

import time
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
                help="untilted stack (black particles in IMAGIC format)")
        parser.add_option("-t",dest="tilted",type="string",metavar="FILE",
                help="tilted stack (black particles in IMAGIC format)")
        parser.add_option("-m",dest="model",type="string",metavar="FILE",
                help="3D model(s) for alignment (Single MRC volume)")
        parser.add_option("-p",dest="param",type="string", metavar="FILE",
                help="Parameter file with refinement info (free_param.par)")
        parser.add_option("-c",dest="ctf1",type="string", metavar="FILE",
                help="Per-particle CTF info for UNTILTED particles")
        parser.add_option("-e",dest="ctf2",type="string", metavar="FILE",
                help="Per-particle CTF info for TILTED particles")
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
	ctf1 = params['ctf1']
	ctf2 = params['ctf2']

	#Get current working directory
	script = sys.argv[0]
	cwd = '%s/lib' %(script[:-26])

	#Get parameter info: angular step
	p = open(param,'r')
	a = 'angular' 
	angl = grep(a,p)
	aL = angl.split()
	ang = aL[2]
	
	#Shift
	p = open(param,'r')
	s = 'shift'
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
	p = open(param,'r')
	r = 'radius'
	radiu = grep(r,p)
	radi = radiu.split()
	rad = radi[2]
	ri = str(float(rad)*float(pix))

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
	incr1 = incr
	incr2 = incr

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
        pp1 = 'psi_step'
        pp2 = grep(pp1,p18)
        pp3 = pp2.split()
        psi = pp3[2]

        p18 = open(param,'r')
        pp1 = 'ncheck'
        pp2 = grep(pp1,p18)
        pp3 = pp2.split()
        ncheck = pp3[2]

		
	#Prepare stack for EMAN2 refinement 
        print '\n'
        print 'Converting stack into correct MRC format'
        print '\n'
        
	if debug is True:
                        print '%s/imagic_to_freeHand2.py -f %s --total=%s --box=%s' %(cwd,untilt,tot,box)

        cmd = '%s/imagic_to_freeHand2.py -f %s --total=%s --box=%s' %(cwd,untilt,tot,box)
        subprocess.Popen(cmd,shell=True).wait()

	#Convert parameter file into freeHand format        
	print '\n'
        print 'Converting CTF info'
        print '\n'

	cmd = '%s/make_f_space_param2.py %s %s' %(cwd,ctf1,mag)
	if debug is True:
		print cmd
	subprocess.Popen(cmd,shell=True).wait()

        #Run refinement
        print '\n'                
	print 'Running search_fspace'                
	print '         ncheck = %s' %(ang)                
	print '         psi_step  = %s' %(sx)                
	print '         RI = %s' %(ri)                
	print '         RMAX1 = %s' %(pix)                
	print '         RMAX2 = %s' %(rad)                
	print '		SH = %s' %(sx)
	print '\n'                

	i = 1

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

                cmd = '%s/../search_refine_fspace/search_fspace_v1_02.csh %s %s %s %s %s %s %s %s %s %s %s %s %s %s.mrc %s_format %s/../search_refine_fspace' %(cwd,first,last,pix,snr,cs,volt,model,ncheck,psi,sx,ri,min_res,max_res,untilt[:-4],ctf1[:-4],cwd)

		if debug is True:
			print cmd
		subprocess.Popen(cmd,shell=True)

		i = i + float(incr)

	#Waiting loop:
	time.sleep(10)
	i = 1

	while i <= int(tot):

		n = str(i + float(incr2) -1)
		n = n[:-2]
		if i == 1:
			first = str(i)
		else:
			first = str(i)
			first = first[:-2]
			incr = incr2
		if float(n) > int(tot):
			incr = int(incr2) - (int(n)-int(tot))
			n = str(tot)

		f = '%s_format_%s_%s' %(ctf1[:-4],first,n)
		fsize = 0
		while fsize == 0:
	
			test = os.path.exists(f)
			if test is False:
				print "Error, %s does not exist; must be bug in refine_fspace" %(f)
				sys.exit()
	
			fsize = os.path.getsize(f)
			time.sleep(1)

			if fsize > 0:
				i = i + float(incr2) 
	tmp = str(float(incr2)-1)

	cmd = '%s/../search_refine_fspace/combine_parfiles.py %s_format_ %s %s' %(cwd,ctf1[:-4],tot,tmp[:-2])
	if debug is True:
		print cmd
	subprocess.Popen(cmd,shell=True).wait()	

	cmd = 'rm %s_format_*_*' %(ctf1[:-4])
	subprocess.Popen(cmd,shell=True).wait()

	print '\n'                
	print 'Converting files into free-hand format'                
	print '\n'                

	cmd = '%s/fspace_param_consolidate.py %s_format_merge.par %s' %(cwd,ctf1[:-4],ctf2)
	subprocess.Popen(cmd,shell=True).wait()
	
	#Convert tilted particles to 3D-MRC format                

	if debug is True:
		print '%s/imagic_to_freeHand2.py -f %s --total=%s --box=%s' %(cwd,tilt,tot,box)
			
	cmd = '%s/imagic_to_freeHand2.py -f %s --total=%s --box=%s' %(cwd,tilt,tot,box)                
	subprocess.Popen(cmd,shell=True).wait()

               
	#Run Free-Hand test                
	info = linecache.getline('%s_format' %(ctf2[:-4]),4)                

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

		cmd = '%s/fastFreeHand_wrapper.csh %s %s %s %s %s.mrc %s %s_format %s %s %s %s %s %s %s %s %s %s %s %s model00 %s' %(cwd,pix,snr,cs,volt,tilt[:-4],model,ctf2[:-4],angSearch,min_res,max_res,str(float(pix)*float(rad)),first,last,incr,mag,df1,df2,astig,str(iteration),calc)

		if debug is True:
			print cmd	

                if debug is False:
                       	subprocess.Popen(cmd,shell=True)

       	        i = i + float(incr)
               	iteration = iteration + 1


if __name__ == "__main__":     
	imagicroot = getIMAGICPath()     
	params=setupParserOptions()     
	main(params)


