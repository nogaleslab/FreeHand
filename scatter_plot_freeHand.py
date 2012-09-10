#!/usr/bin/env python

import optparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#=========================
def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog -p <peaks file> --ang=<float> --lim=<float> --inc=<float>")
        parser.add_option("-p",dest="peaks",type="string",metavar="FILE",
                help="Peaks text file (column 2 = X, column 3 = Y values)")
        parser.add_option("--ang",dest="ang",type="int", metavar="INT",
                help="Expected tilt angle difference")
	parser.add_option("--lim",dest="lim",type="int", metavar="INT",
                help="Angular limit of plot")
	parser.add_option("--inc",dest="inc",type="int", metavar="INT",
                help="Percent of particles included")
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


def scatter(data,lim,tilt,include):

        tiltX = tilt[0]
        tiltY = tilt[1]

        loadData = np.loadtxt(data)
        x = loadData[:,2]
        y = loadData[:,1]

        #Center peak vales at (0,0)
        centx = np.subtract(x,lim)
        centy = np.subtract(y,lim)

        #Calculate distance of peaks from expected angle
        dist = []
        for i in xrange(len(loadData[:,1])):
	        rx = centx[i]
        	ry = centy[i]

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

#Need this at the end for the parse commands
if __name__ == "__main__":
	params=setupParserOptions()
	tiltCenter = [params['ang'],0]
	scatter(params['peaks'],params['lim'],tiltCenter,params['inc'])
