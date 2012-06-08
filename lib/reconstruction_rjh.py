# Adapted Richard J. Hall 11/19/2010 rjhall@berkeley.edu
#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#
from EMAN2_cppwrap import *
from global_def import *

def recons_from_fftvol(size, fftvol, weight, symmetry, npad):
        params = {"size":size, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "weight":weight}
        r = Reconstructors.get("nn4", params)
        r.setup()
        return r.finish(True)


def prepare_recons(data, symmetry, myid, main_node_half, half_start, step, index, finfo=None, npad = 4):
	from random     import randint
	from utilities  import reduce_EMData_to_root
	from mpi        import mpi_barrier, MPI_COMM_WORLD
	nx = data[0].get_xsize()
#	from memorymonitor import MemoryMonitor
	
#	memory_mon = MemoryMonitor('rjhall')

	fftvol_half = EMData()
	weight_half = EMData()
	half_params = {"size":nx, "npad":npad, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = Reconstructors.get( "nn4", half_params )
#	print memory_mon.usage()
	half.setup()
#	print memory_mon.usage()

	group = -1
	for i in xrange(half_start, len(data), step):
		if(index >-1 ):  group = data[i].get_attr('group')
		if(group == index):
			if( data[i].get_attr_default('active',1) == 1):
				xform_proj = data[i].get_attr( "xform.projection" )
				half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	reduce_EMData_to_root(fftvol_half, myid, main_node_half)
	reduce_EMData_to_root(weight_half, myid, main_node_half)
	
	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = randint(0, 1000000)
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi_barrier(MPI_COMM_WORLD)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:  return fftvol_half_file, weight_half_file

	return None, None


def prepare_recons_ctf(nx, data, snr, symmetry, myid, main_node_half, half_start, step, finfo=None, npad = 4):
	from random     import randint
	from utilities  import reduce_EMData_to_root
	from mpi        import mpi_barrier, MPI_COMM_WORLD

        
	fftvol_half = EMData()
	weight_half = EMData()
	half_params = {"size":nx, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = Reconstructors.get( "nn4_ctf", half_params )
	half.setup()

	for i in xrange(half_start, len(data), step):
		if( data[i].get_attr_default('active',1) == 1):
			xform_proj = data[i].get_attr( "xform.projection" )
			half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
        	finfo.write( "begin reduce half\n" )
        	finfo.flush()

	reduce_EMData_to_root(fftvol_half, myid, main_node_half)
	reduce_EMData_to_root(weight_half, myid, main_node_half)
	
	
	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = randint(0, 1000000) 
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi_barrier(MPI_COMM_WORLD)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:
		return fftvol_half_file, weight_half_file

	return None,None

def get_image_size( imgdata, myid ):
	from mpi import mpi_gather, mpi_bcast, MPI_COMM_WORLD, MPI_INT
	nimg = len(imgdata)

        nimgs = mpi_gather( nimg, 1, MPI_INT, 1, MPI_INT, 0, MPI_COMM_WORLD )

        if myid==0:
		src = -1
		for i in xrange( len(nimgs) ):
			if int(nimgs[i]) > 0 :
				src = i
				break
		if src==-1:
			return 0
	else:
		src = -1

	size_src = mpi_bcast( src, 1, MPI_INT, 0, MPI_COMM_WORLD )

	if myid==int(size_src[0]):
		assert nimg > 0
		size = imgdata[0].get_xsize()
	else:
		size = -1

	nx = mpi_bcast( size, 1, MPI_INT, size_src[0], MPI_COMM_WORLD )
	return int(nx[0])


def rec3D_MPI(data, snr, symmetry, mask3D, fsc_curve, myid, main_node = 0, rstep = 1.0, odd_start=0, eve_start=1, finfo=None, index=-1, npad = 4):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept 
          in the memory, computes reconstruction and through odd-even, in order to get the resolution
	'''
	import os
	from statistics import fsc_mask
	from utilities  import model_blank, reduce_EMData_to_root, get_image, send_EMData, recv_EMData
	from random     import randint
	from mpi        import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	
	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005


        if index !=-1 :
		grpdata = []
		for i in xrange( len(data) ):
		    if data[i].get_attr( 'group' ) == index:
		    	    grpdata.append( data[i] )
        	imgdata = grpdata
        else:
		imgdata = data
	nx = get_image_size( imgdata, myid )
	if nx==0:
		ERROR("Warning: no images were given for reconstruction, this usually means there is an empty group, returning empty volume","rec3D",0)
		return model_blank( 2, 2, 2 ), None, model_blank(nx,nx,nx), model_blank(nx,nx,nx)
	
	fftvol_odd_file,weight_odd_file = prepare_recons_ctf(nx, imgdata, snr, symmetry, myid, main_node_odd, odd_start, 2, finfo, npad)
	fftvol_eve_file,weight_eve_file = prepare_recons_ctf(nx, imgdata, snr, symmetry, myid, main_node_eve, eve_start, 2, finfo, npad)
	del imgdata

	if nproc == 1:
		fftvol = get_image(fftvol_odd_file)
		weight = get_image(weight_odd_file)
		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad)

		fftvol = get_image(fftvol_eve_file)
		weight = get_image(weight_eve_file)
		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad)
                
		fscdat = fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)

		fftvol = get_image( fftvol_odd_file )
		fftvol_tmp = get_image(fftvol_eve_file)
		fftvol += fftvol_tmp
		fftvol_tmp = None

		weight = get_image( weight_odd_file )
		weight_tmp = get_image(weight_eve_file)
		weight += weight_tmp
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )
 
		return volall,fscdat,volodd,voleve
  
	if nproc == 2:
		if myid == main_node_odd:
			fftvol = get_image( fftvol_odd_file )
			weight = get_image( weight_odd_file )
			volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad)
			voleve = recv_EMData(main_node_eve, tag_voleve)
			fscdat = fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
		else:
			assert myid == main_node_eve
			fftvol = get_image( fftvol_eve_file )
			weight = get_image( weight_eve_file )
			voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad)
			send_EMData(voleve, main_node_odd, tag_voleve)

		if myid == main_node_odd:
			fftvol = get_image( fftvol_odd_file )
			fftvol_tmp = recv_EMData( main_node_eve, tag_fftvol_eve )
			fftvol += fftvol_tmp
			fftvol_tmp = None

			weight = get_image( weight_odd_file )
			weight_tmp = recv_EMData( main_node_eve, tag_weight_eve )
			weight += weight_tmp
			weight_tmp = None
        	
			volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad)
			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
 
			return volall,fscdat,volodd,voleve
		else:
			assert myid == main_node_eve
			fftvol = get_image( fftvol_eve_file )
			weight = get_image( weight_eve_file )
			send_EMData(fftvol, main_node_odd, tag_fftvol_eve )
			send_EMData(weight, main_node_odd, tag_weight_eve )
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )
			return model_blank(nx,nx,nx), None, model_blank(nx,nx,nx), model_blank(nx,nx,nx)

	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = get_image( fftvol_odd_file )
		send_EMData(fftvol, main_node_eve, tag_fftvol_odd )

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = get_image( weight_odd_file )
		send_EMData(weight, main_node_all, tag_weight_odd )

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad)
		del fftvol, weight
		voleve = recv_EMData(main_node_eve, tag_voleve)
		fscdat = fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
		volall = recv_EMData(main_node_all, tag_volall)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
		return volall,fscdat,volodd,voleve

	if myid == main_node_eve:
		ftmp = recv_EMData(main_node_odd, tag_fftvol_odd)
		fftvol = get_image( fftvol_eve_file )
		Util.add_img( ftmp, fftvol )
		send_EMData(ftmp, main_node_all, tag_fftvol_eve )
		del ftmp

		weight = get_image( weight_eve_file )
		send_EMData(weight, main_node_all, tag_weight_eve )

		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad)
		send_EMData(voleve, main_node_odd, tag_voleve)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );

		return model_blank(nx,nx,nx), None, model_blank(nx,nx,nx), model_blank(nx,nx,nx)


	if myid == main_node_all:
		fftvol = recv_EMData(main_node_eve, tag_fftvol_eve)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = recv_EMData(main_node_odd, tag_weight_odd)
		weight_tmp = recv_EMData(main_node_eve, tag_weight_eve)
		Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad)
		send_EMData(volall, main_node_odd, tag_volall)

		return model_blank(nx,nx,nx),None, model_blank(nx,nx,nx), model_blank(nx,nx,nx)

        return model_blank(nx,nx,nx),None, model_blank(nx,nx,nx), model_blank(nx,nx,nx)

#==========================================
def hsymVols(volodd,voleve,volall,hparams):
	"""
	use volall to find helical parameters
	apply helical parameters to all the volumes
	"""

	hfile = hparams['hfile']
	apix = hparams['apix']
	# find symmetry for full volume (not even or odd)
	lmask = hparams['lmask']
	# save a copy of the unsymmetrized volume
	nosymf = hparams['nosymout']
	volall.write_image(nosymf,-1)
	rot,rise = findHsym(volall,hfile,apix,hparams['isearch'],hparams['osearch'])
	seam = hparams['seam']

	# apply symmetry to all volumes
	volodd = applyHsym(volodd,rot,rise,apix,lmask,seam=seam)
	voleve = applyHsym(voleve,rot,rise,apix,lmask,seam=seam)
	volall = applyHsym(volall,rot,rise,apix,lmask,seam=seam)

	return volodd,voleve,volall

#==========================================
def rec3D_MPI_noCTF(data, symmetry, mask3D, fsc_curve, myid, main_node = 0, rstep = 1.0, odd_start=0, eve_start=1, finfo=None, index = -1, npad = 4, hparams=None):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept in the memory 
	  Computes reconstruction and through odd-even, in order to get the resolution
	  if index > -1, projections should have attribute group set and only those whose group matches index will be used in the reconstruction
	    this is for multireference alignment
	'''
	import os
	from statistics import fsc_mask
	from utilities  import model_blank, reduce_EMData_to_root, get_image,send_EMData, recv_EMData
	from random     import randint
	from mpi        import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	nproc = mpi_comm_size(MPI_COMM_WORLD)
       
	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005
 
        nx = data[0].get_xsize()

        fftvol_odd_file,weight_odd_file = prepare_recons(data, symmetry, myid, main_node_odd, odd_start, 2, index, finfo, npad)
        fftvol_eve_file,weight_eve_file = prepare_recons(data, symmetry, myid, main_node_eve, eve_start, 2, index, finfo, npad) 
	
	if nproc == 1:
		fftvol = get_image( fftvol_odd_file )
		weight = get_image( weight_odd_file )
		volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

		fftvol = get_image( fftvol_eve_file )
		weight = get_image( weight_eve_file )
		voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

		fftvol = get_image( fftvol_odd_file )
		Util.add_img( fftvol, get_image(fftvol_eve_file) )

		weight = get_image( weight_odd_file )
		Util.add_img( weight, get_image(weight_eve_file) )

		volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

		# if helical, find & apply symmetry to volume
		if hparams is not None:
			volodd,voleve,volall = hsymVols(volodd,voleve,volall,hparams)

		fscdat = fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)

		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );
		return volall,fscdat,volodd,voleve

	if nproc == 2:
		if myid == main_node_odd:
			fftvol = get_image( fftvol_odd_file )
			weight = get_image( weight_odd_file )
			volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			voleve = recv_EMData(main_node_eve, tag_voleve)
		else:
			assert myid == main_node_eve
			fftvol = get_image( fftvol_eve_file )
			weight = get_image( weight_eve_file )
			voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			send_EMData(voleve, main_node_odd, tag_voleve)

		if myid == main_node_odd:
			fftvol = get_image( fftvol_odd_file )
			fftvol_tmp = recv_EMData( main_node_eve, tag_fftvol_eve )
			Util.add_img( fftvol, fftvol_tmp )
			fftvol_tmp = None

			weight = get_image( weight_odd_file )
			weight_tmp = recv_EMData( main_node_eve, tag_weight_eve )
			Util.add_img( weight, weight_tmp )
			weight_tmp = None
			volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

			# if helical, find & apply symmetry to volume
			if hparams is not None:
				volodd,voleve,volall = hsymVols(volodd,voleve,volall,hparams)
			fscdat = fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)

			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
			return volall,fscdat,volodd,voleve
		else:
			assert myid == main_node_eve
			fftvol = get_image( fftvol_eve_file )
			send_EMData(fftvol, main_node_odd, tag_fftvol_eve )

			weight = get_image( weight_eve_file )
			send_EMData(weight, main_node_odd, tag_weight_eve )
			import os
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );
			return model_blank(nx,nx,nx), None, model_blank(nx,nx,nx), model_blank(nx,nx,nx)
	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = get_image( fftvol_odd_file )
		send_EMData(fftvol, main_node_eve, tag_fftvol_odd )

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = get_image( weight_odd_file )
		send_EMData(weight, main_node_all, tag_weight_odd )

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		del fftvol, weight
		voleve = recv_EMData(main_node_eve, tag_voleve)
		volall = recv_EMData(main_node_all, tag_volall)

		# if helical, find & apply symmetry to volume
		if hparams is not None:
			volodd,voleve,volall = hsymVols(volodd,voleve,volall,hparams)
		fscdat = fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)

		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
		return volall,fscdat,volodd,voleve

	if myid == main_node_eve:
		ftmp = recv_EMData(main_node_odd, tag_fftvol_odd)
		fftvol = get_image( fftvol_eve_file )
		Util.add_img( ftmp, fftvol )
		send_EMData(ftmp, main_node_all, tag_fftvol_eve )
		del ftmp

		weight = get_image( weight_eve_file )
		send_EMData(weight, main_node_all, tag_weight_eve )

		voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		send_EMData(voleve, main_node_odd, tag_voleve)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );

		return model_blank(nx,nx,nx), None, model_blank(nx,nx,nx), model_blank(nx,nx,nx)


	if myid == main_node_all:
		fftvol = recv_EMData(main_node_eve, tag_fftvol_eve)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = recv_EMData(main_node_odd, tag_weight_odd)
		weight_tmp = recv_EMData(main_node_eve, tag_weight_eve)
		Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		send_EMData(volall, main_node_odd, tag_volall)

		return model_blank(nx,nx,nx),None, model_blank(nx,nx,nx), model_blank(nx,nx,nx)


	return model_blank(nx,nx,nx), None, model_blank(nx,nx,nx), model_blank(nx,nx,nx)

#===========================
def findHsym(vol,hfile,apix,isearch,osearch):
	"""
	use helical search program to find helical symmetry
	"""
	import sys,shutil,subprocess
	from alignment  import helios

	rot,rise = readHsym(hfile)
	orad=osearch/apix
	irad=isearch/apix
	print apix,rise,rot,0.5,orad,irad
	v,newrise,newrot=helios(vol,apix,rise,rot,0.5,orad,irad)

	# store helical parameters
	f=open(hfile,'a')
	f.write("%6.3f  %6.3f\n"%(newrot,newrise))
	f.close()

	return newrot,newrise

#===========================
def applyHsym(vol,rot,rise,apix,lmask,seam=False):
	"""
	apply helical symmetry based on results from Egelman search
	"""
	import shutil,subprocess

	# convert rise to pixels
	nx = vol.get_xsize()
	rise/=apix
	sym=int(round(360.0/abs(rot)))

	if seam is True:
		# apply protofilament symmetry
		sumvol = vol.copy()
		pfoffset=int(sym/2)
		for pnum in range(-pfoffset,sym-pfoffset):
			if pnum==0:
				continue
			ang = rot*pnum
			trans = -(rise*pnum)
			#print pnum, ang, trans
			t = Transform({"type":"spider","psi":ang})
			t.set_trans(0,0,trans)
			volcopy = vol.process("xform",{"transform":t})
			sumvol.add(volcopy)

		# mask in z direction before vertical symmetry
		#sumvol.process_inplace("mask.gaussian.nonuniform",{"radius_x":nx,"radius_y":nx,"radius_z":lmask/apix})	

	else:
		# if applying helical symmetry use himpose
		# hpar file must be in same directory
		shutil.copy(hfile,"tmphpar.spi")
	
		volrot="tmpvol_rt.spi"
		# volume must be rotated for himpose
		t = Transform({"type":"spider","theta":90.0,"psi":90.0})
		volcopy = vol.process("xform",{"transform":t})
		volcopy.write_image(volrot,0,EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
		### EMAN2 SPIDER CONVERSION NOT WORKING PROPERLY
		### WORKAROUND WITH EMAN1
		emancmd = "proc3d %s %s spidersingle"%(volrot,volrot)
		subprocess.Popen(emancmd,shell=True).wait()

		# find himpose exe
		h_exe="himpose"
		himposeexe = subprocess.Popen("which "+h_exe, shell=True, stdout=subprocess.PIPE).stdout.read().strip()
		if not os.path.isfile(himposeexe):
			printError("executable '%s' not found, make sure it's in your path"%h_exe)
		hcmd = "%s %s tmphpar.spi symvol.spi %.3f 2.0 %.1f"%(himposeexe,volrot,apix,int(nx*apix/2)-4)
		print hcmd
		subprocess.Popen(hcmd,shell=True).wait()

		sumvol = EMData("symvol.spi")	
		os.remove(volrot)
		os.remove("tmphpar.spi")
		os.remove("symvol.spi")

	# don't vertically symmetrize 15pf microtubules:
#	if sym!=15:
#		sumvol = verticalCopies(sumvol,sym,rot,rise,seam)
	return sumvol

#===========================
def smart_add(vol1,vol2):
	# add two volumes seamlessly
	nx = vol1.get_xsize()
	sumvol = EMData(nx,nx,nx)
	sumvol.to_zero()
	for i,j,k in ((i,j,k) for i in range(nx) for j in range(nx) for k in range(nx)):
		if vol1[i,j,k] == 0 or vol2[i,j,k] == 0:
			sumvol.set(i,j,k,vol1[i,j,k]+vol2[i,j,k])
		else:
			# at the overlapping area
			sumvol.set(i,j,k,max(vol1[i,j,k],vol2[i,j,k]))
	return sumvol

#===========================
def readHsym(hfile):
	# get helical parameters from file
	f = open(hfile)
	lines = f.readlines()
	f.close()
	pars = lines[-1].strip().split()
	if len(pars) < 2:
		printError("Error: last line of hfile must be helical parameters")
	return float(pars[0]),float(pars[1])

#===========================
def applyHsym_seam(vol,wedgemask,hfile,apix,lmask):
	"""
	apply pseudo-helical symmetry and mask densities on the seam
	"""
	import shutil,subprocess

	# get helical params
	rot,rise=readHsym(hfile)

	# convert rise to pixels
	nx = vol.get_xsize()
	rise/=apix
	sym=int(round(360.0/abs(rot)))

	# apply protofilament symmetry
	vol *= wedgemask
	sumvol = vol.copy()
	pfoffset=int(sym/2)
	for pnum in range(-pfoffset,sym-pfoffset):
		if pnum==0:
			continue
		ang = -(rot*pnum)
		trans = rise*pnum
		#print pnum, ang, trans
		t = Transform({"type":"spider","psi":ang})	
		t.set_trans(0,0,trans)
		volcopy = vol.process("xform",{"transform":t})
		seammaskcopy = wedgemask.process("xform",{"transform":t})
		seammaskcopy.process_inplace("threshold.binary",{"value":0.00001})
		volcopy *= seammaskcopy
		try:
			sumvol.addsmart(volcopy)
		except:
			sumvol = smart_add(sumvol,volcopy)
		
	sumvol.process_inplace("normalize")
	return sumvol
	
#===========================
def verticalCopies(vol,sym,rot,rise,seam):
	rpt = 1.0
	strt = 1
	end = 3
	# if reconstruction has a seam, use dimer as repeat 
	if seam is True:
		rpt = 2.0
		strt = -5
		end = 5
	# apply vertical symmetry
	vertsum = vol.copy()
	for pnum in range (strt,end):
		if pnum==0:
			continue
		smallrot=(sym*abs(rot))%360
		if abs(rot) < (360.0/sym):
			smallrot-=360
		ang = pnum*rpt*(smallrot/3)
		if rot < 0:
			ang*=-1
		trans = -(sym/3.0*rpt)*rise*pnum
		t = Transform({"type":"spider","psi":ang})
		t.set_trans(0,0,trans)
		volcopy = vol.process("xform",{"transform":t})
		vertsum.add(volcopy)
	vertsum.process_inplace("normalize")

	return vertsum

#===========================
def recons_ctf_from_fftvol(size, fftvol, weight, snr, symmetry, npad, weighting=1):
        params = {"size":size, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol, "weight":weight, "weighting":weighting}
        r = Reconstructors.get("nn4_ctf", params)
        r.setup()
        return r.finish(True)

#===========================
def align3Dvols(refvol,vol,apix):
	"""
	Aligns vol to refvol, returns the aligned volume
	Search limited to rotation & translation along z axis,
	and only within 1 protofilament monomer
	"""

	# normalize the vols:
	refvol.process_inplace("normalize")
	vol.process_inplace("normalize")

	rt={}
	rt['maxshift']=40/6.02
	rt['stepdelta']=0
	rt['stepphi']=1
	rt['stepx']=0
	rt['stepy']=0
	rt['stepz']=1

	ali = vol.align("refine.3d",refvol,rt)
	alignvol=vol.process("xform",{"transform":ali["xform.align3d"]})
	
	return alignvol

