#!/bin/csh
#psize,wgh,cs,akv (pixel seize in A, amplitude contrast weight, Cs in mm, voltake in kV, e.g. "2.8,0.07,2.0,200"
#ctfexppart,ctfexpmod (CTF exponent to apply to experimental particles, map projections: e.g. "1,0")
#infile1 (input particle stack, e.g. tilted to -15 degrees "tilted-15.stk")
#infile2 (input reference map, e.g. "model3d.mrc")
#parfile (input parameter file, e.g. for particles tilted to +15 degrees "tilted+15.par")
#outfile1 (output freehand plots for each particle; will need to average them later, e.g. "plots.map")
#rotxstart,rotxstop,rotystart,rotystop (rotation plot limits in degrees, e.g. "-45,45,-45,45")
#rmax1,rmax2,ri (low res freq, high res freq, particle radius in Angstroms, e.g. "300,20,170")
#ifirst,ilast  (first particle and last particle, e.g. "1,100")
#outstyle (p/P=phase residual, c/C=correlation coefficient, e.g "C")
# CTF information for particles in the input stack (infile1) in frealign format.
# NIN,ABSMAGPIN,FILMIN,DFMID1IN,DFMID2IN,ANGASTIN,IMORE
#
# NIN = number of particle in (e.g. "50")
# ABSMAGPIN = Absolute magification of particles (e.g. "50000")
# FIMLIN = film number corresponding to NIN particles (e.g. "1234")
# DFMID1IN = Defocus midpoint 1 in Angstroms (e.g. "35000")
# DFMID2IN = Defocus midpoint 2 in Angstroms (e.g. "36000")
# ANGASTIN = Angle of astigmatism for CTF in degrees (e.g. "35.00")
# IMORE = another line of CTF information to follow ("1"=yes, "2"=no)
#

set pix = $1
set amp = $2
set cs = $3
set volt = $4
set stack = $5
set model = $6
set param = $7
set ang = $8
set min = $9
set max = $10
set rad = $11
set first = $12
set last = $13
set num = $14
set mag = $15
set df1 = $16
set df2 = $17
set astig = $18
set iter = $19
set name = $20
set cc = $21

time ~michael/BATCHLIB/freeHand/fastfreehand_v1_01.exe << eot
$pix,$amp,$cs,$volt
1,0
$stack
$model
$param
${name}_plots_CC_v101_${iter}.mrc
-${ang},$ang,-${ang},$ang
$min,$max,$rad
$first,$last
$cc
$num,$mag,1,$df1,$df2,$astig,0
eot

