#!/usr/bin/env python

import glob
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
        parser.set_usage("%prog -p <parameter file> -m <model> -o <output folder>")
        parser.add_option("-p",dest="param",type="string", metavar="FILE",
                help="Parameter file with refinement info (free_param.par)")
        parser.add_option("-m",dest="model",type="string", metavar="FILE",
                help="Model(s) file used for free-hand test")
	parser.add_option("-o",dest="out",type="string", metavar="FILE",
                help="Output folder name")
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

def peak(stack,tot,cent):
	spifile = "currentSpiderScript.spi"
       	if os.path.isfile(spifile):
        	os.remove(spifile)
       	spi=open(spifile,'w')
	spicmd="SD IC NEW\n"
	spicmd+="incore\n"
	spicmd+="2,%s\n" %(tot)
	spicmd+="do lb1 [part] = 1, %s\n" %(tot)
	spicmd+="\n"
	spicmd+="PK [x] [y]\n"
	spicmd+="%s@{******[part]}\n"%(stack[:-4])
	spicmd+="(1,0)\n"
	spicmd+="[newX] = %s +[x]\n" %(cent)
	spicmd+="[newY] = %s +[y]\n" %(cent)
	spicmd+="SD IC [part] [newX] [newY]\n"
	spicmd+="incore\n"
	spicmd+="lb1\n"
	spicmd+="SD IC COPY\n"
	spicmd+="incore\n"
	spicmd+="%s_peak\n" %(stack[:-4])
	spicmd+="SD ICE\n"
	spicmd+="incore\n"
	spicmd+="en d\n"
	runSpider(spicmd)

def runSpider(lines):
       spifile = "currentSpiderScript.spi"
       if os.path.isfile(spifile):
               os.remove(spifile)
       spi=open(spifile,'w')
       spi.write("MD\n")
       spi.write("TR OFF\n")
       spi.write("MD\n")
       spi.write("VB OFF\n")
       spi.write("MD\n")
       spi.write("SET MP\n")
       spi.write("(0)\n")
       spi.write("\n")
       spi.write(lines)

       spi.write("\nEN D\n")
       spi.close()
       spicmd = "spider spi @currentSpiderScript"
       spiout = subprocess.Popen(spicmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr.read()
       output = spiout.strip().split()
       if "ERROR" in output:
               print "Spider Error, check 'currentSpiderScript.spi'\n"
               sys.exit()
       # clean up
       os.remove(spifile)
       if os.path.isfile("LOG.spi"):
               os.remove("LOG.spi")
       resultf = glob.glob("results.spi.*")
       if resultf:
               for f in resultf:
                       os.remove(f)

def getIMAGICPath():
        ### get the imagicroot directory
        impath = subprocess.Popen("env | grep IMAGIC_ROOT", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        imagicpath = impath.replace("IMAGIC_ROOT=","")
        if imagicpath != '/opt/qb3/imagic-110326':
                        print "imagic/110326 was not found, make sure it is in your path"
                        sys.exit()

def getCCP4Path():        
        ### get the openmpi directory        
        ccp4path = subprocess.Popen("env | grep CCP4_PATH", shell=True, stdout=subprocess.PIPE).stdout.read().strip()        

        if ccp4path:                
                ccp4path = ccp4path.replace("CCP4_PATH=","")      
        if os.path.exists(ccp4path):                        
                return ccp4path        
        print "ccp4 is not loaded, make sure it is in your path"        
        sys.exit()

def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            return match.string

def getEMANPath():        
        ### get the imagicroot directory        
        emanpath = subprocess.Popen("env | grep EMAN2DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()        

        if emanpath:                
                emanpath = emanpath.replace("EMAN2DIR=","")                
        if os.path.exists(emanpath):                        
                return emanpath        
        print "EMAN2 was not found, make sure eman2/2.05 is in your path"        
        sys.exit()

def main2(params):
	
	param = params['param']
	out = params['out']
	debug = params['debug']
	model = params['model']
        #Free hand angular search
        p8 = open(param,'r')
        fs1 = 'freeHand_ang_search'
        fs2 = grep(fs1,p8)
        fs3 = fs2.split()
        angSearch = fs3[2]

        p18 = open(param,'r')
        fs1 = 'calc'
        fs2 = grep(fs1,p8)
        fs3 = fs2.split()
        calc = fs3[2]

        #Free hand angular search
        p9 = open(param,'r')
        fs1 = 'num_mod'
        fs2 = grep(fs1,p9)
        fs3 = fs2.split()
        mods = float(fs3[2])

	#Merge stacks:
	m = 0

	script = sys.argv[0]
	cwd = '%s/lib' %(script[:-22])
	
	while m < mods:

		num = len(glob.glob('model%02d_plots*.mrc'%(m)))

		i = 1
		while i <= int(num):

			cmd = '%s/mrc_to_im.b model%02d_plots_CC_v101_%s.mrc' %(cwd,m,i)
			if debug is True:
				print cmd
			subprocess.Popen(cmd,shell=True).wait()

			cmd = 'proc2d model%02d_plots_CC_v101_%s.img model%02d_plots_CC_v101_merge.img' %(m,i,m)
			if debug is True:
				print cmd
			subprocess.Popen(cmd,shell=True).wait()

			i = i + 1

		cmd = '%s/im_to_mrc.b model%02d_plots_CC_v101_merge.img' %(cwd,m)
		subprocess.Popen(cmd,shell=True).wait()

	        cmd = 'cp %s/totsumstack.exe .' %(cwd)
	        subprocess.Popen(cmd,shell=True).wait()

		cmd = 'cp %s/totsumstack_mult.csh .' %(cwd)
                subprocess.Popen(cmd,shell=True).wait()
	
	        cmd = './totsumstack_mult.csh model%02d' %(m) 
	        subprocess.Popen(cmd,shell=True).wait()

	        cmd = 'rm totsumstack.exe totsumstack_mult.csh'
	        subprocess.Popen(cmd,shell=True).wait()
	
		if calc is 'C':

			line1 = (float(angSearch)*2)/5
			line = line1/2
			if debug is True:
				print '%s = float(%s*2)/5' %(line1,angSearch)
				print '%s = %s / 2' %(line,line1)
				print '%s/npo_CC_wrap_mult.csh %s model%02d %s' %(cwd,str(float(angSearch)*2),m,line)

		        cmd = '%s/npo_CC_wrap_mult.csh %s model%02d %s' %(cwd,str(float(angSearch)*2),m,line)
		        subprocess.Popen(cmd,shell=True).wait()

			cmd = '%s/mrc_to_spi.b model%02d_plots_CC_v101_merge.mrc' %(cwd,m)
			subprocess.Popen(cmd,shell=True).wait()
	
			cmd = 'mkdir %s' %(out)
	       		subprocess.Popen(cmd,shell=True).wait()
	
			cmd = 'mv model%02d_plots_CC_v101* %s/' %(m,out)
		        subprocess.Popen(cmd,shell=True).wait()

			cmd = 'mv model%02d_frehand_CC.ps %s/' %(m,out)
		        subprocess.Popen(cmd,shell=True).wait()

	       		cmd = 'mv model%02d_averageplot_CC_v101.mrc %s/' %(m,out)
		        subprocess.Popen(cmd,shell=True).wait()

			cmd = 'cp %s %s/' %(model,out)
			subprocess.Popen(cmd,shell=True).wait()
	
		        cmd = 'cp %s %s/' %(param,out)
		        subprocess.Popen(cmd,shell=True).wait()
	
			tot = EMUtil.get_image_count('%s/model%02d_plots_CC_v101_merge.img' %(out,m)) 	
			n = int(angSearch)+1
			stack = '%s/model%02d_plots_CC_v101_merge.spi' %(out,m)
			peak(stack,tot,n)
			m = m + 1

                if calc is 'P':

                        line1 = (float(angSearch)*2)/5
                        line = line1/2
                        if debug is True:
                                print '%s = float(%s*2)/5' %(line1,angSearch)
                                print '%s = %s / 2' %(line,line1)
                                print '%s/npo_CC_wrap_mult_phase.csh %s model%02d %s' %(cwd,str(float(angSearch)*2),m,line)

                        cmd = '%s/npo_CC_wrap_mult_phase.csh %s model%02d %s' %(cwd,str(float(angSearch)*2),m,line)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = '%s/mrc_to_spi.b model%02d_plots_CC_v101_merge.mrc' %(cwd,m)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'mkdir %s' %(out)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'mv model%02d_plots_CC_v101* %s/' %(m,out)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'mv model%02d_frehand_CC.ps %s/' %(m,out)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'mv model%02d_averageplot_CC_v101.mrc %s/' %(m,out)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'cp %s %s/' %(model,out)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'cp %s %s/' %(param,out)
                        subprocess.Popen(cmd,shell=True).wait()

                        tot = EMUtil.get_image_count('%s/model%02d_plots_CC_v101_merge.img' %(out,m))
                        n = int(angSearch)+1
                        stack = '%s/model%02d_plots_CC_v101_merge.spi' %(out,m)
                        peak(stack,tot,n)
                        m = m + 1

		cmd = 'mv *00.* %s' %(out)	
		subprocess.Popen(cmd,shell=True).wait()

		cmd = 'rm -r logfile* test.img test.hed model??*.mrc refine_eman2 z.plot start.hdf *_prep.img *_prep.hed '
  		subprocess.Popen(cmd,shell=True).wait()

		cmd = "cp %s/find_peaks_freeHand.spi %s" %(cwd,out)
		subprocess.Popen(cmd,shell=True).wait()

if __name__ == "__main__":
        imagicroot = getIMAGICPath()
        getCCP4Path()
	getEMANPath()
	from EMAN2 import *
        params=setupParserOptions()
        main2(params)

