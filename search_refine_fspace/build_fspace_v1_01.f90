!-------------------------------------------------------
! build_fspace Version 1.01 Dated 060212
!-------------------------------------------------------
! JLR 08/08 Original code
! JLR 01/09 Made addimage a subroutine
! JLR 04/09 Some optimizations added, revised to take multiple stack/parameter files
!           Norming volume converted to DP.  Workingvol overwritten by normalized volume
!           Fixed bug in FSC
! JLR 08/10 Check memory allocations for maps as they are done. 
!           Fixed bug in renormalize map subroutine (define norm=0 at beginning)
!           Fixed bug in FSC (ony fill curve out to Nyquist instead of overflow)
!           Allow > and < thresholding
! JLR 04/11 New modular design implemented, completed 11/11
! JLR 02/12 Fixed output titles on map

INCLUDE 'fspace_subsandfuncs_vI.f90'

PROGRAM Build_fspace
USE Consts
USE Imageinfo
USE CTFmod
USE Fouriervolume, ONLY: Addimage, Renormalizemap, Fsctest, Moveorigin2d

IMPLICIT NONE

! Image

INTEGER                       :: part,nparts,xpos,ypos
INTEGER                       :: ifirst,ilast
REAL                          :: dmin,dmax,dmean,radius
INTEGER                       :: nxyz(3),nxyzst(3),mxyz(3),mode
INTEGER                       :: vox,sec,lin
INTEGER                       :: nin, ifilmin,imore,ilist
REAL                          :: amagp,presa,cell(6)
REAL                          :: dfmid1last,dfmid2last,angastlast
REAL                          :: absmag,dfmid1,dfmid2,angast,angasti
INTEGER                       :: film
CHARACTER(100)                :: line

INTEGER                       :: i,j,k

CHARACTER(80)                 :: outfile1,outfile2,title,format1,format2,maskfile
CHARACTER(80), ALLOCATABLE    :: instack(:),parfile(:)
CHARACTER(24)                 :: dat

REAL                          :: norm
INTEGER                       :: normterms
REAL, ALLOCATABLE             :: image(:,:), map(:,:,:), mask(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE   :: partc(:,:),imagec(:,:)
DOUBLE COMPLEX, ALLOCATABLE   :: workingvol(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE   :: workingvolA(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE   :: workingvolB(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE   :: tempvolcA(:,:,:),tempvolcB(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: tempvolA(:,:,:),tempvolB(:,:,:)
INTEGER                       :: memoryerror


DOUBLE PRECISION, ALLOCATABLE :: normingvol(:,:,:),normingvolA(:,:,:),normingvolB(:,:,:)


DOUBLE PRECISION, ALLOCATABLE :: ctfx(:,:),ctfx_wrap(:,:)
REAL, ALLOCATABLE             :: fsc(:)
INTEGER                       :: usemask

LOGICAL                :: newctf
REAL                   :: psize,wgh,cs,akv,wgh1,wgh2,wl,wc1,wc2,thetatr
!REAL                   :: ctf
REAL                   :: ctfexppart
INTEGER                :: now(3) ! Time information

! Interpolation
!REAL                   :: sinclut(2000)
REAL                   :: pshftr
REAL                   :: kx,ky,kz,kxplane,kyplane,kzplane,arg(3)
INTEGER                :: kxs,kxf,kys,kyf,kzs,kzf,kxvol,kyvol,kzvol,&
                          kxarray,kyarray,kzarray,kxwrap,kywrap,kzwrap
INTEGER                :: ksum 

! Orientation parameters
REAL                   :: insert(6)
REAL                   :: psi,theta,phi,shx,shy
INTEGER                :: shifts(2)
REAL                   :: psiout,thetaout,phiout
REAL                   :: cpsi,spsi,cthe,sthe,cphi,sphi
REAL                   :: rmax1,rmax2,rmax1sq,rmax2sq,radiussq

! Stack management
INTEGER                :: stack,stackparts,stackposition,nstacks,described

! Thresholding
INTEGER                :: included, excluded
REAL                   :: thresh
CHARACTER(1)           :: thresharg

! FFTW arrays
DOUBLE PRECISION, ALLOCATABLE :: image_fftw(:,:)
DOUBLE COMPLEX, ALLOCATABLE   :: imagec_fftw(:,:)
DOUBLE PRECISION, ALLOCATABLE :: vol_fftw(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE   :: volc_fftw(:,:,:)
REAL*8                        :: plan,plan2d,plan3df,plan3dr

! FFTW parameters
INCLUDE '/usr/include/fftw3.f'


!!$!-------------------------------------------------
!!$! Calculate a sinc LUT
!!$!-------------------------------------------------
!!$WRITE (*,*) "calculating SINC lookup table"
!!$DO j=1,2000
!!$  sinclut(j)=SIN(FLOAT(j)*PI/180.0)/(FLOAT(j)*PI/180.0)
!!$END DO
!!!CALL Calculatesinclut

!-------------------------------------------------
! get file information and open files
!-------------------------------------------------
READ (5,*)  psize,wgh,cs,akv
READ (5,*)  ctfexppart
READ (5,*)  rmax2
READ (5,*)  ifirst,ilast,nstacks
READ (5,*)  outfile1,usemask,maskfile,thresharg,thresh ! map,FSCmask

ALLOCATE(instack(nstacks))
ALLOCATE(parfile(nstacks))

DO i=1,nstacks
   READ (5,'(A80)')  instack(i) ! particle stack
   READ (5,'(A80)')  parfile(i) ! input Parameters
ENDDO

! Stream 1 = particle stack input
! Stream 2 = map output
! Stream 3 = parameter file input
! Stream 4 = FSC mask
! Stream 5 = Reserved for executing file
! Stream 6 = Reserved for screen output
! Stream 7 = Output FSC file

! Open 1st particle stack
CALL Imopen(1,instack(1),"RO")
CALL Irdhdr(1,nxyz,mxyz,mode,dmin,dmax,dmean)
nx     = nxyz(1)
ny     = nxyz(2)
nz     = nxyz(1)
stackparts = nxyz(3)
nparts = (ilast-ifirst)+1
CALL GenOptimizationConstants
! Restate resolution limit in terms of pixel radii
rmax2=(psize*nx)/rmax2
rmax2sq=rmax2**2
rmax1sq=0


!-----------------------------------------
! Open map
!-----------------------------------------
CALL Imopen(2,outfile1,"new")

!-------------------------------------------------
! Allocate arrays
!-------------------------------------------------
radius = 0.5 * nx
! allocate matrices for input and output maps
ALLOCATE (image(nx,ny))                   ; image=0
ALLOCATE (map(nx,ny,nz))                  ; map=0
ALLOCATE (mask(nx,ny,nz))                 ; mask=0
ALLOCATE (partc(nx/2+1,ny))               ; partc=0
ALLOCATE (ctfx(nx/2+1,ny))                ; ctfx=0
ALLOCATE (ctfx_wrap(nx/2+1,ny))           ; ctfx_wrap=0
ALLOCATE (image_fftw(nx,ny))              ; image_fftw=0
ALLOCATE (imagec_fftw(nx/2+1,ny))         ; imagec_fftw=0
ALLOCATE (imagec(nx/2+1,ny))              ; imagec=0
ALLOCATE (vol_fftw(nx,ny,nz))             ; vol_fftw=0
ALLOCATE (volc_fftw(nx/2+1,ny,nz))        ; volc_fftw=0
ALLOCATE (fsc(nx/2+1))                    ; fsc=0
! Large memory allocations
ALLOCATE (workingvol(nx/2+1,ny,nz),stat=memoryerror)       
IF (memoryerror/=0) STOP 'memory allocation error'; workingvol=0
ALLOCATE (normingvol(nx/2+1,ny,nz),stat=memoryerror)
IF (memoryerror/=0) STOP 'memory allocation error'; normingvol=0
ALLOCATE (workingvolA(nx/2+1,ny,nz),stat=memoryerror)      
IF (memoryerror/=0) STOP 'memory allocation error'; workingvolA=0
ALLOCATE (normingvolA(nx/2+1,ny,nz),stat=memoryerror)      
IF (memoryerror/=0) STOP 'memory allocation error'; normingvolA=0
ALLOCATE (workingvolB(nx/2+1,ny,nz),stat=memoryerror)      
IF (memoryerror/=0) STOP 'memory allocation error'; workingvolB=0
ALLOCATE (normingvolB(nx/2+1,ny,nz),stat=memoryerror)      
IF (memoryerror/=0) STOP 'memory allocation error'; normingvolB=0
ALLOCATE (tempvolcA(nx/2+1,ny,nz),stat=memoryerror)        
IF (memoryerror/=0) STOP 'memory allocation error'; tempvolcA=0
ALLOCATE (tempvolcB(nx/2+1,ny,nz),stat=memoryerror)        
IF (memoryerror/=0) STOP 'memory allocation error'; tempvolcB=0
ALLOCATE (tempvolA(nx,ny,nz),stat=memoryerror)             
IF (memoryerror/=0) STOP 'memory allocation error'; tempvolA=0
ALLOCATE (tempvolB(nx,ny,nz),stat=memoryerror)             
IF (memoryerror/=0) STOP 'memory allocation error'; tempvolB=0
! allocate parameterfile variable


!-----------------------------------------
! Load in FSC Maskfile
!-----------------------------------------
IF (usemask==1) THEN
  CALL Imopen(4,maskfile,"RO")
  CALL Irdhdr(4,nxyz,mxyz,mode,dmin,dmax,dmean)
  DO sec=1,nz
    DO lin=1,ny
      CALL Imposn(4,sec-1,lin-1) 
      CALL Irdlin (4,mask(:,lin,sec),*999)
    END DO 
  END DO
  CALL Imclose(4)
ELSE IF (usemask==0) THEN
  mask=1
ELSE
  WRITE (*,*) "The number that preceeds the name of the maskfile should be '1' if the"
  WRITE (*,*) "maskfile is to be used, or '0' if no mask is to be used for FSC calculation."
  STOP
END IF
  

!-----------------------------------------
! Get CTF Information 
!-----------------------------------------
format1='(I7,5F8.3,F8.0,I6,2F9.1,F8.2,F7.2)'
format2='(I7,5F8.3,F8.0,I6,2F9.1,F8.2,F7.2,I3.2,A1,I2.2,A1,I2.2)'



!-------------------------------------------------
! Get to correct starting point for stack and parameter file
!-------------------------------------------------

!stackparts = nxyz(3)
! determined earlier
described = stackparts
stack=1
stackposition=ifirst
DO
  IF (stackposition<=stackparts) THEN  
    EXIT
  ELSE
    CALL Imclose(1)
    stack=stack+1
    stackposition=stackposition-stackparts
    CALL Imopen(1,instack(stack),"RO")
    CALL Irdhdr(1,nxyz,mxyz,mode,dmin,dmax,dmean)
    stackparts = nxyz(3)
    described=described+stackparts
  END IF
END DO

! Open the correct parameter file
OPEN(unit=3,file=parfile(stack),status="unknown")
! Read past comment lines
DO
  READ (3,'(A100)') line
  IF (line(1:23).eq."C           PSI   THETA") EXIT
END DO

!-------------------------------------------------
! Main Loop over all particles
!-------------------------------------------------
AKV=1000.0*AKV
WL=12.3/SQRT(AKV+AKV**2/(10.0**6.0))
wgh1=SQRT(1.0-wgh**2)
wgh2=wgh
thetatr=wl/(psize*nx)   

CALL Dfftw_plan_dft_r2c_2d(plan2d,nx,ny,image_fftw,imagec_fftw,fftw_measure)
dfmid1last=0.000001
dfmid2last=0.000001
angastlast=0.000001

stackposition=stackposition-1
included=0; excluded=0
Eachpart: DO part=ifirst,ilast
  stackposition=stackposition+1
  !-------------------------------------------------
  ! Open a new stack and parameter file if necessary
  !-------------------------------------------------
  !!
  IF (part>described) THEN
    stack=stack+1
    CALL Imclose(1)
    CALL Imopen(1,instack(stack),"RO")
    CALL Irdhdr(1,nxyz,mxyz,mode,dmin,dmax,dmean)
    stackparts = nxyz(3)
    stackposition=1
    described=described+stackparts
    CLOSE(3)
    OPEN(unit=3,file=parfile(stack),status="unknown")
    DO
      READ (3,'(A100)') line
      IF (line(1:23).eq."C           PSI   THETA") EXIT
    END DO
  END IF
  !-------------------------------------------------
  ! put image values into array parts(:,:,:)
  !-------------------------------------------------
  !write (*,*) "about to put particles into array"
  DO lin=1,ny
    CALL Imposn(1,stackposition-1,lin-1) 
    CALL Irdlin (1,image(:,lin),*999)    
  END DO
  !-------------------------------------------------
  ! read in the relevant line of the parameter file
  !-------------------------------------------------
  DO
    READ(3,format1) ilist,psi,theta,phi,shx,shy,&
                    amagp,film,dfmid1,dfmid2,angasti,presa
    IF (ilist==stackposition) EXIT
  END DO
  angast=angasti*pi/180
  psi=psi*pi/180
  theta=theta*pi/180
  phi=phi*pi/180
  !!shx=2*pi*shx/nx !! Changed fspace_subsandfuncs to accept pixel shifts for shx and shy
  !!shy=2*pi*shy/ny
  !-------------------------------------------------
  ! Generate a new CTF multiplier if necessary (not wraparound order)
  !-------------------------------------------------
  IF (dfmid1last.ne.dfmid1.or.dfmid2last.ne.dfmid2.or.angastlast.ne.angasti) THEN
    WRITE (*,'(A19,I8)') "New CTF at particle",part
    dfmid1last=dfmid1
    dfmid2last=dfmid2
    angastlast=angasti
    ctfx=0
    CALL GenCTFmultiplier(ctfx,ctfx_wrap,nxby2,nyby2,ny,rmax1sq,rmax2sq,cs,wl,wgh1,wgh2,dfmid1,dfmid2,angast,thetatr)
  END IF
  !!WRITE (*,*) "Working on particle:",part
  !-------------------------------------------------
  ! Obtain FT of image and apply CTF
  !-------------------------------------------------
  IF ((thresharg=='<'.and.presa<thresh).or.(thresharg=='>'.and.presa>thresh)) THEN
    included=included+1
    image_fftw(:,:)=DBLE(image(:,:))  
    CALL Dfftw_execute(plan2d)
    CALL Moveorigin2d(imagec_fftw,nxby2,nxby2+1,ny)
    ! Take image out of wraparound order
    imagec(1:nxby2+1,1:nyby2-1)=imagec_fftw(1:nxby2+1,nyby2+2:ny)
    imagec(1:nxby2+1,nyby2:ny)=imagec_fftw(1:nxby2+1,1:nyby2+1)
    IF (IAND(part,1).eq.0) THEN
      CALL Addimage(workingvolA,normingvolA,imagec,nx,ny,nz,phi,theta,psi,&
           shx,shy,rmax2sq,ctfx,ctfexppart)
    ELSE
      CALL Addimage(workingvolB,normingvolB,imagec,nx,ny,nz,phi,theta,psi,&
           shx,shy,rmax2sq,ctfx,ctfexppart)
    END IF
  ELSE
    excluded=excluded+1
  END IF

  CALL Itime(now)
END DO Eachpart


!-------------------------------------------------
! Prepare and transform Fourier volume
!-------------------------------------------------
CALL Dfftw_plan_dft_c2r_3d(plan3dr,nx,ny,nz,volc_fftw,vol_fftw,fftw_measure)
CALL Dfftw_plan_dft_r2c_3d(plan3df,nx,ny,nz,vol_fftw,volc_fftw,fftw_measure)


! Combine the two half maps
workingvol=workingvolA+workingvolB
normingvol=normingvolA+normingvolB
CALL Renormalizemap(workingvolA,normingvolA,nx,ny,nz)
CALL Renormalizemap(workingvolB,normingvolB,nx,ny,nz)


!-------------------------------------------------
! Apply realspace masks to working volumes
!-------------------------------------------------
WRITE(*,*) "Moving half maps into wrap-around order for reverse FT"
! Move map into wraparound order and move RS origin to centre for masking
DO kzarray=1,nz
   IF(kzarray.lt.nzby2) kzwrap=kzarray+nzby2+1
   IF(kzarray.ge.nzby2) kzwrap=kzarray-nzby2+1
   DO kyarray=1,ny      
      IF(kyarray.lt.nyby2) kywrap=kyarray+nyby2+1
      IF(kyarray.ge.nyby2) kywrap=kyarray-nyby2+1
      DO kxarray=1,nxby2+1
         ksum=(kxarray-1)+(kyarray-1)+(kzarray-1)
         pshftr=1.0
         IF (MOD(ksum,2).ne.0) pshftr=-1.0
         kxwrap=kxarray
         tempvolcA(kxwrap,kywrap,kzwrap)=workingvolA(kxarray,kyarray,kzarray)*pshftr
         tempvolcB(kxwrap,kywrap,kzwrap)=workingvolB(kxarray,kyarray,kzarray)*pshftr
      END DO
   END DO
END DO

! Apply mask if needed
IF (usemask==1) THEN
  volc_fftw=tempvolcA
  CALL Dfftw_execute(plan3dr)
  tempvolA=vol_fftw*nxnynzinv
  WRITE(*,*) "Reverse FT'd first half map"
  volc_fftw=tempvolcB
  CALL Dfftw_execute(plan3dr)
  tempvolB=vol_fftw*nxnynzinv
  WRITE(*,*) "Reverse FT'd second half map"
   vol_fftw=tempvolA*DBLE(mask)
  CALL Dfftw_execute(plan3df)
  tempvolcA=volc_fftw
  vol_fftw=tempvolB*DBLE(mask)
  CALL Dfftw_execute(plan3df)
  tempvolcB=volc_fftw
END IF

! Put maps back into wraparound order
DO kzarray=1,nz
   IF (kzarray.lt.nzby2) kzwrap=kzarray+nzby2+1
   IF (kzarray.ge.nzby2) kzwrap=kzarray-nzby2+1
   DO kyarray=1,ny
      IF (kyarray.lt.nyby2) kywrap=kyarray+nyby2+1
      IF (kyarray.ge.nyby2) kywrap=kyarray-nyby2+1
      DO kxarray=1,nxby2+1
         kxwrap=kxarray
         workingvolA(kxarray,kyarray,kzarray)=tempvolcA(kxwrap,kywrap,kzwrap)
         workingvolB(kxarray,kyarray,kzarray)=tempvolcB(kxwrap,kywrap,kzwrap)
      END DO
   END DO
END DO

! perform Fourier Shell Correlation test
CALL fsctest(workingvolA,workingvolB,fsc,nx,ny,nz)
! Output Fourier Shell Correlation

OPEN(unit=7,file='fsc.txt',status="unknown")
WRITE (*,'(A6,A7,A8)') "C Ring","Res","FSC"
WRITE (7,'(A6,A7,A8)') "C Ring","Res","FSC"
DO kxarray=2,nx/2+1
  IF (kxarray>rmax2) EXIT
  WRITE (*,'(I7,F7.2,F8.3)') kxarray,psize*nx/(kxarray-1),fsc(kxarray)
  WRITE (7,'(I7,F7.2,F8.3)') kxarray,psize*nx/(kxarray-1),fsc(kxarray)
END DO

IF (usemask==0) THEN
  WRITE (*,'(A53)') "C No masking of models for FSC calculation was performed"
  WRITE (7,'(A56)') "C No masking of models for FSC calculation was performed"
ELSE
  WRITE (*,'(A53)') "C Masking of models for FSC calculation was performed"
  WRITE (7,'(A53)') "C Masking of models for FSC calculation was performed"
END IF

WRITE (*,'(A19,I10,A21,I10)') "Particles included:",included,"  Particles excluded:",excluded

CALL Renormalizemap(workingvol,normingvol,nx,ny,nz)

! Move map into wraparound order and move origin to centre
DO kxarray=1,nx/2+1
  DO kyarray=1,ny
    DO kzarray=1,nz
      ksum=(kxarray-1)+(kyarray-1)+(kzarray-1)
      pshftr=1.0
      IF (MOD(ksum,2).ne.0) pshftr=-1.0
      kxwrap=kxarray
      IF (kyarray<=ny/2+1) kywrap=kyarray+ny/2-1 ! wrap around in y   
      IF (kyarray>ny/2+1) kywrap=kyarray-ny/2-1 ! positive y
      IF (kzarray<=nz/2+1)   kzwrap=kzarray+nz/2-1 ! wrap around in z
      IF (kzarray>nz/2+1) kzwrap=kzarray-nz/2-1 ! positive z
      volc_fftw(kxarray,kyarray,kzarray)=workingvol(kxwrap,kywrap,kzwrap)*pshftr
    END DO
  END DO
END DO

CALL Dfftw_execute(plan3dr)
DO vox=1,nx
  DO lin=1,ny
    DO sec=1,nz
      !map(vox,lin,sec)=REAL(vol_fftw(vox,lin,sec))*mask(vox,lin,sec)/(nx*ny*nz) ! Build with mask
      map(vox,lin,sec)=REAL(vol_fftw(vox,lin,sec))*nxnynzinv ! Build w/o mask
    END DO
  END DO
END DO

CALL Dfftw_destroy_plan(plan2d)
CALL Dfftw_destroy_plan(plan3dr)
CALL Dfftw_destroy_plan(plan3df)

!-------------------------------------------------
! write out map
!-------------------------------------------------

WRITE (*,*) "Writing output map"
dmin = MINVAL(map)
dmax = MAXVAL(map)
dmean = SUM(map)*nxnynzinv

title="3D map created by build_fspace"
nxyz=(/nx,ny,nz/)
nxyzst=(/0,0,0/)
mxyz=(/nx,ny,nz/)
cell=(/nx*psize,ny*psize,nz*psize,90.0,90.0,90.0/)
WRITE (*,*) "Writing output map"
CALL Itrhdr(2,1)
CALL Ialsiz(2,nxyz,nxyzst)
CALL Ialsam(2,mxyz)
CALL Ialcel(2,cell)
CALL Iwrhdr(2,title,0,dmin,dmax,dmean)
DO sec=1,nz
  DO lin=1,ny
    CALL Imposn(2,sec-1,lin-1) 
    CALL Iwrlin(2,map(:,lin,sec))
  END DO
END DO

CALL Imclose(1)
CALL Imclose(2)
CLOSE(3)


!-------------------------------------------------
! end
!-------------------------------------------------


STOP "Normal termination of Build_fspace"
999 STOP "End of File Read Error"
END PROGRAM Build_fspace
