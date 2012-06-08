!-------------------------------------------------------
! Fastfrehand Version 1.01 Dated 310112
!-------------------------------------------------------
! JLR 10/08
! JLR 05/09
! JLR 03/11
! JLR 07/11 New modular design implemented
! JLR 01/12 Bugfix: optimization constants generated before nz defined; code cleanup
INCLUDE 'fspace_subsandfuncs_modules_v1_01.f90'
PROGRAM Fastfrehand
USE Consts
USE Imageinfo
USE CTFmod
USE Fouriervolume, ONLY: Extractimage_noshifts, Moveorigin2d
USE OtherFoms, ONLY: Pres
IMPLICIT NONE

! Counters, constants, positions
INTEGER                           :: i,j,k
INTEGER                           :: part,nparts,xpos,ypos
INTEGER                           :: ifirst,ilast
REAL                              :: radius
INTEGER                           :: vox,sec,lin
INTEGER                           :: nin, ifilmin, imore

! Physical/CTF parameters
REAL                              :: absmagpin,dfmid1in,dfmid2in,angastin,ccbestin
REAL                              :: dfmid1last,dfmid2last,angastlast
REAL, ALLOCATABLE                 :: absmag(:),dfmid1(:),dfmid2(:),angast(:),ccbest(:)
INTEGER, ALLOCATABLE              :: film(:)
REAL                              :: psize,wgh,cs,akv,wgh1,wgh2,wl,wc1,wc2,thetatr
REAL                              :: ctfexppart,ctfexpmod
INTEGER                           :: ilist
REAL                              :: amagp

! Image parameters
REAL                              :: dmin,dmax,dmean
INTEGER                           :: nxyz(3),mxyz(3),cell(6),mode

! File names and titles
CHARACTER(80)                     :: infile1,infile2,parfile,outfile1,outfile2,title,format1,format2
CHARACTER(100)                    :: line
CHARACTER(1)                      :: outstyle

! Image and map arrays
REAL, ALLOCATABLE                 :: parts(:,:,:), projs(:,:,:),map(:,:,:)
DOUBLE PRECISION, ALLOCATABLE     :: rimask(:,:),ctfx(:,:),ctfx_wrap(:,:)
LOGICAL, ALLOCATABLE              :: ccfmask(:,:)
DOUBLE COMPLEX, ALLOCATABLE       :: partsc(:,:,:), imagec(:,:)
DOUBLE COMPLEX, ALLOCATABLE       :: workingvol(:,:,:),workingslice(:,:)

! Map wrapping/unwrapping
DOUBLE PRECISION                  :: pshftr
INTEGER                           :: kxarray,kyarray,kzarray,kxwrap,kywrap,kzwrap
INTEGER                           :: ksum


! Cross-correlation 
REAL                          :: CCC,ssref
REAL, ALLOCATABLE             :: sumsquaredpart(:)
INTEGER                       :: ncheckx,nchecky
INTEGER                       :: rotxcheck, rotycheck,rotxstart,rotxstop,rotystart,rotystop
INTEGER                       :: plotxpos,plotypos
REAL                          :: taxis

! Rotation matrices
REAL*8                        :: r1(3,3),r2(3,3),rg(3,3),r12(3,3),r2in(3,3),r1in(3,3),re(3,3),rp(3,3),rg2(3,3)
REAL*8                        :: va(3),vb(3),vc(3),vd(3),ve(3),vf(3),axisi(3),axisg(3),cp(3)
REAL*8                        :: vg(3),vh(3),vi(3),vj(3)
REAL                          :: DET
REAL*8                        :: rotx,roty,rotz,tang,xcomp,ycomp

! Search parameters
REAL*8                        :: psicheck,thetacheck,phicheck
REAL                          :: extract(6)
REAL, ALLOCATABLE             :: psi(:),theta(:),phi(:),shx(:),shy(:)
INTEGER                       :: shifts(2)
REAL                          :: shxbest,shybest
REAL                          :: cpsi,spsi,cthe,sthe,cphi,sphi
REAL                          :: rmax1,rmax2,rmax1sq,rmax2sq,ri

! FFTW arrays
DOUBLE PRECISION, ALLOCATABLE :: image_fftw(:,:)
DOUBLE COMPLEX, ALLOCATABLE   :: imagec_fftw(:,:)
DOUBLE PRECISION, ALLOCATABLE :: vol_fftw(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE   :: volc_fftw(:,:,:)
REAL*8                        :: plan2df,plan2dr,plan3d

! Optimization variables
DOUBLE PRECISION              :: ccnorm

INCLUDE '/usr/include/fftw3.f'
!-------------------------------------------------
! get file information and open files
!-------------------------------------------------
READ (5,*)        psize,wgh,cs,akv
READ (5,*)        ctfexppart,ctfexpmod
READ (5,'(A80)')  infile1  ! input particle stack
READ (5,'(A80)')  infile2  ! input map
READ (5,'(A80)')  parfile  ! input parameter file
READ (5,'(A80)')  outfile1 ! output frehand plots
READ (5,*)        rotxstart,rotxstop,rotystart,rotystop
READ (5,*)        rmax1,rmax2,ri
READ (5,*)        ifirst,ilast
READ (5,*)        outstyle

! Datastream 
! 1=particle stack
! 2=map
! 3=parameter file
! 4=output frehandplots

! Open particle stack
WRITE (*,*) "Opening particle stack",ifirst,ilast
CALL Imopen(1,infile1,"RO")
CALL Irdhdr(1,nxyz,mxyz,mode,dmin,dmax,dmean)
nx      = nxyz(1)
ny      = nxyz(2)
nparts  = (ilast-ifirst)+1
ncheckx = rotxstop-rotxstart+1
nchecky = rotystop-rotystart+1 

! Restate resolution limits in terms of pixel radii
rmax1=(psize*nx)/rmax1
IF (rmax1<1) rmax1=1.01 ! Never use the origin of the FT for alignment
rmax2=(psize*nx)/rmax2
ri=ri/psize

rmax1sq=rmax1**2
rmax2sq=rmax2**2

!-----------------------------------------
! Open map
!-----------------------------------------

WRITE (*,*) "Getting reference model information"
CALL Imopen(2,infile2,"RO")
CALL Irdhdr(2,nxyz,mxyz,mode,dmin,dmax,dmean)
nz     = nxyz(3)
IF (nxyz(1).ne.nx) STOP 'File sizes do not match'

! nx, ny, and nz all defined
CALL GenOptimizationConstants
!-------------------------------------------------
! Allocate arrays
!-------------------------------------------------

WRITE (*,*) "Allocating arrays"
radius = 0.5 * nx
! allocate matrices for input and output maps
ALLOCATE (map(nx,ny,nx))                ; map=0
ALLOCATE (parts(nx,ny,nparts))          ; parts=0
ALLOCATE (sumsquaredpart(nparts))       ; sumsquaredpart=0
ALLOCATE (projs(ncheckx,nchecky,nparts+1)); projs=0
ALLOCATE (ccfmask(nx,ny))               ; ccfmask=0
ALLOCATE (rimask(nx,ny))                ; rimask=0
ALLOCATE (partsc(nxby2plus1,ny,nparts)) ; partsc=0
ALLOCATE (ctfx(nx/2+1,ny))              ; ctfx=0
ALLOCATE (ctfx_wrap(nx/2+1,ny))         ; ctfx_wrap=0
ALLOCATE (absmag(nparts))               ; absmag=0
ALLOCATE (workingvol(nxby2plus1,ny,nz)) ; workingvol=0
ALLOCATE (workingslice(nxby2plus1,ny))  ; workingslice=0
ALLOCATE (image_fftw(nx,ny))            ; image_fftw=0
ALLOCATE (imagec_fftw(nxby2plus1,ny))   ; imagec_fftw=0
ALLOCATE (imagec(nxby2plus1,ny))        ; imagec=0
ALLOCATE (vol_fftw(nx,ny,nz))           ; vol_fftw=0
ALLOCATE (volc_fftw(nxby2plus1,ny,nz))  ; volc_fftw=0
! Parameter file info
ALLOCATE (psi(nparts))               ; psi=0
ALLOCATE (theta(nparts))             ; theta=0
ALLOCATE (phi(nparts))               ; phi=0
ALLOCATE (shx(nparts))               ; shx=0
ALLOCATE (shy(nparts))               ; shy=0
ALLOCATE (ccbest(nparts))            ; ccbest=0
ALLOCATE (dfmid1(nparts))            ; dfmid1=0
ALLOCATE (dfmid2(nparts))            ; dfmid2=0
ALLOCATE (angast(nparts))            ; angast=0
ALLOCATE (film(nparts))              ; film=0

!-------------------------------------------------
! Plan 2d Fourier transforms
!-------------------------------------------------
!WRITE (*,*) "Preparing FFTW Fourier Transforms"
CALL Dfftw_plan_dft_r2c_2d(plan2df,nx,ny,image_fftw,imagec_fftw,fftw_measure)
CALL Dfftw_plan_dft_c2r_2d(plan2dr,nx,ny,imagec_fftw,image_fftw,fftw_measure)
!WRITE (*,*) "Preparing input and output formats"
format1='(I7,5F8.3,F8.0,I6,2F9.1,F8.2,F7.2)'
format2='(I7,5F8.3,F8.0,I6,2F9.1,F8.2,F7.2,I3.2,A1,I2.2,A1,I2.2)'

!-------------------------------------------------
! Prepare a mask for the CCF
! Prepare a mask for the real space image
!-------------------------------------------------
WRITE (*,*) "Preparing CCF and RI masks"
DO lin=1,ny
  DO vox=1,nx
    IF (lin<=ny/2+1) ypos = lin-1      
    IF (lin>ny/2+1)  ypos = lin-ny-1
    IF (vox<=nx/2)   xpos = vox-1      
    IF (vox>nx/2)    xpos = vox-nx-1
    radius=(xpos**2 + ypos**2)**0.5 ! Use REAL and SQRT statements instead
    IF (radius<=5) THEN
      ccfmask(vox,lin)=.true.
    ELSE
      ccfmask(vox,lin)=.false.
    END IF
  END DO
END DO

DO lin=1,ny
  DO vox=1,nx
    xpos=vox-nx/2
    ypos=lin-ny/2
    radius=(xpos**2 + ypos**2)**0.5
    IF (radius<=ri) THEN
      rimask(vox,lin)=1
    ELSE
      rimask(vox,lin)=0
    END IF
  END DO
END DO
!-----------------------------------------
! Read parameter file
!-----------------------------------------
WRITE (*,*) "Reading input parameter file"
OPEN(unit=3,file=parfile,status="old")
DO
  READ (3,'(A100)') line
  IF (line(1:23).eq."C           PSI   THETA") EXIT
END DO
! Skip irrelevant records in parameter file
i=1
DO 
  READ(3,format1) ilist,psi(i),theta(i),phi(i),shx(i),shy(i),amagp,film(i),dfmid1in,dfmid2in,angastin,ccbestin
  IF (ilist<ifirst) THEN
    i=1
  ELSE
   psi(i)=psi(i)
   theta(i)=theta(i)
   phi(i)=phi(i)
   shx(i)=2*pi*shx(i)/nx
   shy(i)=2*pi*shy(i)/ny
   i=i+1
  END IF
  IF (ilist==ilast) EXIT
END DO


!-----------------------------------------
! Get CTF Information 
!-----------------------------------------
WRITE (*,*) "Getting CTF information for all particles"

part=0
DO
  READ(5,*)NIN,ABSMAGPIN,IFILMIN,DFMID1IN,DFMID2IN,ANGASTIN,IMORE
  DO i=1,NIN
    part=part+1
    IF (part>=nparts) EXIT
    film(part)=ifilmin
    dfmid1(part)=dfmid1in
    dfmid2(part)=dfmid2in
    angast(part)=angastin*pi/180
    absmag(part)=absmagpin
    IF (part>=nparts) EXIT
  END DO
  IF (part>=nparts) EXIT
  IF(imore==0) EXIT
END DO

!-------------------------------------------------
! put image values into array parts(:,:,:)
!-------------------------------------------------
WRITE (*,*) "Putting particle images into array"
part=0
DO sec=ifirst,ilast
  part=part+1
  CALL Imposn(1,sec-1,0) 
  CALL Irdsec (1,parts(:,:,part),*999)    
END DO 

WRITE (*,*) "Putting map into array"
DO sec=1,nz
  CALL Imposn(2,sec-1,0) 
  CALL Irdsec (2,map(:,:,sec),*999)
END DO

!-------------------------------------------------
! Transform particle images, apply CTF
!-------------------------------------------------
WRITE (*,*) "Applying CTF to particle images"
AKV=1000.0*AKV
WL=12.3/SQRT(AKV+AKV**2/(10.0**6.0))
wgh1=SQRT(1.0-wgh**2)
wgh2=wgh
thetatr=wl/(psize*nx)  
dfmid1last=0.000001
dfmid2last=0.000001
angastlast=0.000001

DO part=1,nparts
  !-------------------------------------------------
  ! Generate a new CTF multiplier if necessary
  !-------------------------------------------------
  IF (dfmid1last.ne.dfmid1(part).or.dfmid2last.ne.dfmid2(part).or.angastlast.ne.angast(part)) THEN
    dfmid1last=dfmid1(part)
    dfmid2last=dfmid2(part)
    angastlast=angast(part)
    CALL GenCTFmultiplier(ctfx,ctfx_wrap,nxby2,nyby2,ny,rmax1sq,rmax2sq,CS,WL,WGH1,WGH2,DFMID1(part),&
                          DFMID2(part),ANGAST(part),THETATR)
  END IF
  image_fftw(:,:)=DBLE(parts(:,:,part))
  CALL Dfftw_execute(plan2df)
  DO kyarray=1,ny
    DO kxarray=1,nx/2+1
      partsc(kxarray,kyarray,part)=imagec_fftw(kxarray,kyarray)*ctfx_wrap(kxarray,kyarray)**ctfexppart      
    END DO
  END DO
END DO

!-------------------------------------------------
! Get sumsquaredpart(:) for normalized CCF calculation later
! Put CTF applied particle back into parts(:,:,:) and apply mask
!-------------------------------------------------
WRITE (*,*) "Applying RI mask and getting image statistics"
DO part=1,nparts
  imagec_fftw(:,:)=partsc(:,:,part)
  CALL Dfftw_execute(plan2dr)
  image_fftw=image_fftw*rimask/(nx*ny)
  parts(:,:,part)=REAL(image_fftw)
  CALL Dfftw_execute(plan2df) 
  partsc(:,:,part)=imagec_fftw
  DO lin=1,ny
    DO vox=1,ny
      sumsquaredpart(part)=sumsquaredpart(part)+(parts(vox,lin,part))**2
    END DO
  END DO
END DO
!-------------------------------------------------
! Fourier transform map
!-------------------------------------------------
WRITE (*,*) "Preparing 3D Fourier volume"
CALL Dfftw_plan_dft_r2c_3d(plan3d,nx,ny,nz,vol_fftw,volc_fftw,fftw_measure)
vol_fftw=DBLE(map)
CALL Dfftw_execute(plan3d)
!!CALL Dfftw_destroy_plan(plan3d)


!-------------------------------------------------
! Place map FT in working volume (i.e. origin at centre)
!-------------------------------------------------
DO kxarray=1,nxby2plus1
  DO kyarray=1,ny
    DO kzarray=1,nz
      ksum=(kxarray+1)+(kyarray-1)+(kzarray-1)
      pshftr=1.0
      IF (MOD(ksum,2).ne.0) pshftr=-1.0
      kxwrap=kxarray
      IF (kyarray<nyby2)    kywrap=kyarray+nyby2plus1 
      IF (kyarray>=nyby2)   kywrap=kyarray-nyby2+1 
      IF (kzarray<nzby2)    kzwrap=kzarray+nzby2+1 
      IF (kzarray>=nzby2)   kzwrap=kzarray-nzby2+1 
      workingvol(kxarray,kyarray,kzarray)=volc_fftw(kxwrap,kywrap,kzwrap)*pshftr
    END DO
  END DO
END DO


!-------------------------------------------------
! Main Loop over all particles
!-------------------------------------------------
dfmid1last=0.000001
dfmid2last=0.000001
angastlast=0.000001

WRITE (*,*) "Main loop"
Eachpart: DO part=1,nparts
  ! Get a rotation matrix for the previously assigned Euler angles
  r1=0
  CALL EUL2ROT(R1,DBLE(psi(part)),DBLE(theta(part)),DBLE(phi(part)))
  ! check determinant is near 1.0
  DET=0
  DET = (r1(1,1)*r1(2,2)-r1(1,2)*r1(2,1))*r1(3,3) + & 
        (r1(2,1)*r1(3,2)-r1(2,2)*r1(3,1))*r1(1,3) + &
        (r1(3,1)*r1(1,2)-r1(3,2)*r1(1,1))*r1(2,3)
  IF(ABS(DET-1.0).GT.0.0001) STOP 'Determinant not near 1.0 - euler angles 1'
  ! Calculate predicted euler angles for second image
  WRITE (*,*) "Working on particle",part
  !-------------------------------------------------
  ! Generate a new CTF**2 multiplier if necessary
  !-------------------------------------------------
  IF (dfmid1last.ne.dfmid1(part).or.dfmid2last.ne.dfmid2(part).or.angastlast.ne.angast(part)) THEN
    dfmid1last=dfmid1(part)
    dfmid2last=dfmid2(part)
    angastlast=angast(part)
    CALL GenCTFmultiplier(ctfx,ctfx_wrap,nxby2,nyby2,ny,rmax1sq,rmax2sq,CS,WL,WGH1,WGH2,DFMID1(part),&
                          DFMID2(part),ANGAST(part),THETATR)
    WRITE (*,*) "Generating new CTF multiplier for calculated projections"
  ELSE
  END IF
  !-------------------------------------------------
  ! Cycle through all Euler angles
  !-------------------------------------------------
  plotypos=0
  Eachroty: DO rotycheck=rotystart,rotystop
    plotypos=plotypos+1
    plotxpos=0
    WRITE (*,*) "Working on particle",part
    Eachrotx: DO rotxcheck=rotxstart,rotxstop
      plotxpos=plotxpos+1
      rotx=DBLE(rotxcheck)*pi/180
      roty=DBLE(rotycheck)*pi/180
      tang=SQRT(roty**2+rotx**2)
      taxis=ATAN2(rotx,roty) 
      ! this order is requred to match the angplot convention for rotx, roty
      ! a large rotx value (as defined for angplot) results in a large rotation about the y axis
      xcomp=COS(taxis)
      ycomp=SIN(taxis)
      !
      Rg(1,1) = COS(tang)+(xcomp**2)*(1-COS(tang))
      Rg(1,2) = xcomp*ycomp*(1-COS(tang))
      Rg(1,3) = ycomp*SIN(tang)
      !
      Rg(2,1) = xcomp*ycomp*(1-COS(tang))
      Rg(2,2) = COS(tang)+(ycomp**2)*(1-COS(tang))
      Rg(2,3) = -1*xcomp*SIN(tang)
      !
      Rg(3,1) = -1*ycomp*SIN(tang)
      Rg(3,2) = xcomp*SIN(tang)
      Rg(3,3) = COS(tang)
      ! check determinant is near 1.0
      DET=0
      DET = (rg(1,1)*rg(2,2)-rg(1,2)*rg(2,1))*rg(3,3) + & 
            (rg(2,1)*rg(3,2)-rg(2,2)*rg(3,1))*rg(1,3) + &
            (rg(3,1)*rg(1,2)-rg(3,2)*rg(1,1))*rg(2,3)
      IF(ABS(DET-1.0).GT.0.0001) STOP 'Determinant not near 1.0 - goniometer rotation'
      ! Calculate predicted euler angles for second image
      Rp=matmul(R1,Rg)
      ! check determinant is near 1.0
      DET=0
      DET = (rp(1,1)*rp(2,2)-rp(1,2)*rp(2,1))*rp(3,3) + & 
            (rp(2,1)*rp(3,2)-rp(2,2)*rp(3,1))*rp(1,3) + &
            (rp(3,1)*rp(1,2)-rp(3,2)*rp(1,1))*rp(2,3)
      IF(ABS(DET-1.0).GT.0.0001) STOP 'Determinant not near 1.0 - particle 2 predicted alignment'
      psicheck=0;thetacheck=0;phicheck=0
      CALL ROT2EUL(Rp,psicheck,thetacheck,phicheck)
      CALL Extractimage_noshifts(workingvol,workingslice,nx,ny,nz,REAL(phicheck*pi/180),REAL(thetacheck*pi/180),&
                                 REAL(psicheck*pi/180),rmax1sq,rmax2sq)
      !-------------------------------------------------
      ! Apply the CTF multiplier to the extracted slice
      !-------------------------------------------------
      DO lin=1,ny
        DO vox=1,nx/2+1
          workingslice(vox,lin)=workingslice(vox,lin)*(ctfx(vox,lin)**ctfexpmod)
        END DO
      END DO
      ! Unload working slice into wraparound order
      imagec(1:nxby2plus1,nyby2plus2:ny)=workingslice(1:nxby2plus1,1:nyby2minus1)
      imagec(1:nxby2plus1,1:nyby2plus1)=workingslice(1:nxby2plus1,nyby2:ny)
      ! Move origin
      CALL Moveorigin2d(imagec,nxby2,nxby2plus1,ny)
      ssref=0
      ! FT Into realspace
      imagec_fftw=imagec
      CALL Dfftw_execute(plan2dr)
      image_fftw=image_fftw*nxnyinv
      ! Calculate sum-squared pixel value of reference
      ssref=0
      DO lin=1,ny
        DO vox=1,nx
          ssref=ssref+REAL(image_fftw(vox,lin))**2
        END DO
      END DO
      !-----------------------------------------------------------------
      ! Calculate normalized CCF
      !-----------------------------------------------------------------
      DO k=1,ny
        DO j=1,nxby2plus1
          imagec_fftw(j,k)=imagec(j,k)*CONJG(partsc(j,k,part))
        END DO
      END DO
      CALL Dfftw_execute(plan2dr)
      IF (outstyle=='c'.or.outstyle=="C") THEN
        ccnorm=1/DBLE(SQRT(ssref*sumsquaredpart(part)))
        CCC=MAXVAL(image_fftw,ccfmask)*nxnyinv*ccnorm
        projs(plotxpos,plotypos,part)=CCC
      ELSE IF (outstyle=='p'.or.outstyle=="P") THEN
        shifts=MAXLOC(image_fftw,ccfmask)
        shxbest=REAL(shifts(1))
        shybest=REAL(shifts(2))
        IF (shybest<nyby2plus1)  shybest=-1*shybest+1
        IF (shybest>nyby2plus1) shybest=ny-shybest+1
        IF (shxbest<nxby2plus1) shxbest=-1*shxbest+1
        IF (shxbest>nxby2plus1)  shxbest=nx-shxbest+1
        shxbest=-shxbest
        shybest=-shybest
        projs(plotxpos,plotypos,part)=pres(imagec,partsc(:,:,part),shxbest,shybest,nx,ny,rmax1sq,rmax2sq,nxby2,nyby2)*180/pi
      ELSE
        STOP 'Please select output format as Phase Residual or Correlation coefficient ("P" or "C")'
      END IF
    END DO Eachrotx
  END DO Eachroty
  !-------------------------------------------------
  ! Output best orientation for each particle
  !------------------------------------------------- 
END DO Eachpart

!-------------------------------------------------
! write out plots
!-------------------------------------------------

WRITE (*,*) "Writing output frehand plot(s)"
dmin = MINVAL(projs)
dmax = MAXVAL(projs)
dmean = SUM(projs)/(ncheckx*nchecky*nparts)
title="  "
nxyz=(/ncheckx,nchecky,nparts/)
mxyz=(/ncheckx,nchecky,nparts/)
cell=(/ncheckx,nchecky,nparts,90,90,90/)
CALL Imopen(4,outfile1,"NEW")
CALL Icrhdr(4,nxyz,mxyz,2,"   ",1)
CALL Ialsiz(4,nxyz,nxyz) !!!!changed last parameter from nxyzi
CALL Itrcel(4,2)
CALL Ialcel(4,cell)
CALL Irtcel(4,cell)
CALL Iwrhdr(4,title,-1,dmin,dmax,dmean)
DO sec=1,nparts+1
  DO lin=1,nchecky
    CALL Imposn(4,sec-1,lin-1) 
    CALL Iwrlin(4,projs(:,lin,sec))
  END DO
END DO

CALL Imclose(1)
CALL Imclose(2)
CALL Imclose(4)
CLOSE(4)
CALL Dfftw_destroy_plan(plan2df)
CALL Dfftw_destroy_plan(plan2dr)
!-------------------------------------------------
! end
!-------------------------------------------------
STOP "Normal termination of Fastfrehand"
999 STOP "End of File Read Error"
END PROGRAM Fastfrehand


!******************************************************************************
SUBROUTINE EUL2ROT(RT,PSI,THETA,PHI)
!******************************************************************************
!
IMPLICIT NONE
REAL*8 PI
PARAMETER (PI=3.141592654)
!
INTEGER NANG,J,IAXIS,IANGL,NSYM,JSYM
REAL*8 PSI,THETA,PHI,RT(3,3),DM(9)
REAL*8 CPHI,SPHI,CTHE,STHE,CPSI,SPSI
!
CPHI=COS(PHI*PI/180.)
SPHI=SIN(PHI*PI/180.)
CTHE=COS(THETA*PI/180.)
STHE=SIN(THETA*PI/180.)
CPSI=COS(PSI*PI/180.)
SPSI=SIN(PSI*PI/180.)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Below is the standard Eulerian rotation ZYZ
!    The model or its transform is rotated first by PHI around Z, then 
!    by THETA about the new Y, and thirdly by PSI about the new Z.
!
!    The rotation matrix used is  R=R(psi)*R(theta)*R(phi) as in Spider
!
!                 c  s  0          c  0 -s          c  s  0
!                -s  c  0     *    0  1  0     *   -s  c  0
!                 0  0  1          s  0  c          0  0  1
!
!                 about Z          about Y          about Z
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DM(1)=(CPHI*CTHE*CPSI-SPHI*SPSI)
DM(2)=(SPHI*CTHE*CPSI+CPHI*SPSI)
DM(3)=-STHE*CPSI
DM(4)=(-CPHI*CTHE*SPSI-SPHI*CPSI)
DM(5)=(-SPHI*CTHE*SPSI+CPHI*CPSI)
DM(6)=STHE*SPSI
DM(7)=STHE*CPHI
DM(8)=STHE*SPHI
DM(9)=CTHE

RT(1,1)=DM(1)
RT(2,1)=DM(2)
RT(3,1)=DM(3)
RT(1,2)=DM(4)
RT(2,2)=DM(5)
RT(3,2)=DM(6)
RT(1,3)=DM(7)
RT(2,3)=DM(8)
RT(3,3)=DM(9)

RETURN
END

!*****************************************************************************
SUBROUTINE ROT2EUL(RT,PSI,THETA,PHI)
!*****************************************************************************
!  converts a standard rotation matrix into Eulerian angles 
!  similar to CCP4 MATROT subroutine but with Eulerian ZYZ convention.
!
IMPLICIT NONE
REAL*8 PI
PARAMETER (PI=3.141592654)
!
REAL*8 PSI,THETA,PHI,RT(3,3),RT2(3,3),DRT
REAL*8 DET,CPHI,SPHI,CTHE,STHE,CPSI,SPSI,DRAD
!
DRAD=PI/180.0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Below is the standard Eulerian rotation ZYZ
!    The model or its transform is rotated first by PHI around Z, then 
!    by THETA about the new Y, and thirdly by PSI about the new Z.
!
!    The rotation matrix used is  R=R(psi)*R(theta)*R(phi) as in Spider
!
!                 c  s  0          c  0 -s          c  s  0
!                -s  c  0     *    0  1  0     *   -s  c  0
!                 0  0  1          s  0  c          0  0  1
!
!                 about Z          about Y          about Z
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Eulerian angles can be forced without restriction of generality
!    to be located in the interval 0<=THETA<=180 (this is a consequence
!    of the identity operation
!   PSI, THETA, PHI -> PSI+180, THETA, PHI+180 in Eulerian angle space)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  first check determinant is near 1.0
DET = (RT(1,1)*RT(2,2)-RT(1,2)*RT(2,1))*RT(3,3) + & 
      (RT(2,1)*RT(3,2)-RT(2,2)*RT(3,1))*RT(1,3) + &
      (RT(3,1)*RT(1,2)-RT(3,2)*RT(1,1))*RT(2,3)

IF(ABS(DET-1.0).GT.0.0001) STOP 'Determinant not near 1.0'

!  first get cos(theta),theta and sin(theta)
CTHE=MAX(-1.0,MIN(1.0,RT(3,3)))
THETA=ACOS(CTHE)/DRAD
STHE=SQRT(MAX(0.0,1.0-RT(3,3)**2))

!  for theta not equal to 0 or 180, PHI, PSI are unique
IF(ABS(ABS(CTHE)-1.0).GT.0.00000001) THEN
  CPHI=RT(1,3)/STHE
  SPHI=RT(2,3)/STHE
  CPHI=MAX(-1.0,MIN(1.0,CPHI))
  PHI=ACOS(CPHI)
  IF(SPHI.LT.0.0)PHI=2.0*PI-PHI
  PHI=PHI/DRAD
  CPSI=-RT(3,1)/STHE
  SPSI=RT(3,2)/STHE
  CPSI=MAX(-1.0,MIN(1.0,CPSI))
  PSI=ACOS(CPSI)
  IF(SPSI.LT.0.0)PSI=2.0*PI-PSI
  PSI=PSI/DRAD
ELSE
!  for THETA=0/180, PHI and PSI can have an infinite number of values, only
!   [PSI-PHI] is defined, so PHI can be set to zero without restriction
  PHI=0.0
  CPSI=RT(1,1)
  SPSI=RT(2,1)
  CPSI=MAX(-1.0,MIN(1.0,CPSI))
  PSI=ACOS(CPSI)
  IF(SPSI.LT.0.0)PSI=2.0*PI-PSI
  PSI=PSI/DRAD
END IF

!      DRT=0.1
! now check for consistency with EUL2ROT
!      CALL EUL2ROT(RT2,PSI,THETA,PHI)
!      IF(ABS(RT(1,1)-RT2(1,1)).GT.DRT.OR.
!     .	 ABS(RT(1,2)-RT2(1,2)).GT.DRT.OR.
!     .	 ABS(RT(1,3)-RT2(1,3)).GT.DRT.OR.
!     .	 ABS(RT(2,1)-RT2(2,1)).GT.DRT.OR.
!     .	 ABS(RT(2,2)-RT2(2,2)).GT.DRT.OR.
!     .	 ABS(RT(2,3)-RT2(2,3)).GT.DRT.OR.
!     .	 ABS(RT(3,1)-RT2(3,1)).GT.DRT.OR.
!     .	 ABS(RT(3,2)-RT2(3,2)).GT.DRT.OR.
!     .	 ABS(RT(3,3)-RT2(3,3)).GT.DRT) THEN
!      	WRITE(*,*) ' STOP - inconsistent PSI,THETA,PHI =',PSI,THETA,PHI
!      	WRITE(*,20) RT,RT2
!20	FORMAT(' Input rotation matrix '/3(3F12.6/),
!     .	       ' Output rotation matrix'/3(3F12.6/))
!      	STOP
!      ENDIF
!
RETURN
END
