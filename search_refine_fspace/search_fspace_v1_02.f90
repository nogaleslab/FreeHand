!-------------------------------------------------------
! search_fspace Version 1.01 Dated 290612
!-------------------------------------------------------
! Exhaustive projection matching
!-------------------------------------------------------
! JLR 1/09
! JLR 7/11  New modular design implemented
! JLR 11/11 Upgrade completed
! JLR 03/12 Changed to output standard Frealign format parameter file
!           and use an input Frealign parameter file for CTF info
!           Cross correlation coefficient now expressed out of 100
! JLR 06/12 Removed output of orientation information (orientlist)
!           
!-------------------------------------------------------

INCLUDE 'fspace_subsandfuncs_vI.f90'

PROGRAM Search_fspace
USE Consts
USE Imageinfo
USE CTFmod
USE Fouriervolume, ONLY: Extractimage_noshifts, Moveorigin2d

IMPLICIT NONE
! 3D volume and basic parameters !
INTEGER                           :: nparts,nstacks,nprojs
INTEGER                           :: ifirst,ilast
REAL, ALLOCATABLE                 :: map(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE       :: workingvol(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE       :: workingslice(:,:),workingslice2(:,:),workingslice3(:,:)

! File name informations
CHARACTER(80)                     :: inmap,title
CHARACTER(80),ALLOCATABLE         :: instack(:),parfileout(:),parfilein(:)
CHARACTER(80)                     :: format1,format2,format3,format4,format5,format6

! File parameters
REAL                              :: dmin,dmax,dmean,radius,cell(6)
INTEGER                           :: nxyz(3),nxyzst(3),mxyz(3),mode

! Image and FTs of image arrays
REAL, ALLOCATABLE                 :: image(:,:)
DOUBLE PRECISION, ALLOCATABLE     :: particle(:,:)
DOUBLE COMPLEX, ALLOCATABLE       :: particlec(:,:), imagec(:,:)

! FFTW arrays
DOUBLE PRECISION, ALLOCATABLE     :: image_fftw(:,:)
DOUBLE COMPLEX, ALLOCATABLE       :: imagec_fftw(:,:)
DOUBLE PRECISION, ALLOCATABLE     :: vol_fftw(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE       :: volc_fftw(:,:,:)
REAL*8                            :: plan,plan2df,plan2dr,plan3df,plan3dr

! CTF and other mask variables
REAL                              :: absmag,dfmid1,dfmid2,angast,angasti
LOGICAL                           :: newctf
REAL                              :: psize,wgh,cs,akv,wgh1,wgh2,wl,wc1,wc2,thetatr
REAL                              :: ctfexppart,ctfexpmod
DOUBLE PRECISION, ALLOCATABLE     :: rimask(:,:),ctfx(:,:),ctfx_wrap(:,:)
LOGICAL, ALLOCATABLE              :: ccfmask(:,:)
INTEGER                           :: nin,imore
INTEGER                           :: film,ilist
REAL                              :: amagp
REAL                              :: dfmid1last,dfmid2last,angastlast

! Cross-correlation 
REAL                              :: CCC,CCbest,psibest,thetabest,phibest,ssref,sumsquaredpart
REAL                              :: shxbest,shybest
REAL, ALLOCATABLE                 :: sumsquaredref(:)
INTEGER                           :: orientationID
 
! Time variables
INTEGER                           :: now(3)

! Interpolation
REAL                              :: pshftr
REAL                              :: rmax1,rmax2,rmax1sq,rmax2sq

! Array markers
REAL                              :: kx,ky,kz,kxplane,kyplane,kzplane
INTEGER                           :: kxvol,kyvol,kzvol
INTEGER                           :: kxarray,kyarray,kzarray,kxwrap,kywrap,kzwrap
INTEGER                           :: ksum 
INTEGER                           :: i,j,k
INTEGER                           :: vox,sec,lin
INTEGER                           :: part,xpos,ypos

! Stack managements
INTEGER                           :: stack, stackposition, described, stackparts

! Orientation parameters
REAL                              :: h,psi,theta,phi
REAL                              :: lastphi
REAL                              :: psi_step,ri,shiftmax
REAL                              :: psiout,thetaout,phiout
INTEGER                           :: shifts(2),shx,shy
INTEGER                           :: ncheck,quadrant,hemisphere
CHARACTER(100)                    :: line
REAL                              :: jnkr

! Variables to speed up code
DOUBLE PRECISION                  :: ccnorm          !1/(nx*ny),

INCLUDE '/usr/include/fftw3.f'
!-------------------------------------------------
! get startup information 
!-------------------------------------------------

READ (5,*)        psize,wgh,cs,akv
READ (5,*)        ctfexppart,ctfexpmod
READ (5,'(A80)')  inmap      ! input map
READ (5,*)        ncheck,psi_step,shiftmax,ri
READ (5,*)        rmax1,rmax2
READ (5,*)        ifirst,ilast,nstacks

ALLOCATE (instack(nstacks))
ALLOCATE (parfilein(nstacks))
ALLOCATE (parfileout(nstacks))

DO i=1,nstacks
  READ (5,'(A80)')  instack(i)  ! input particle stack
  READ (5,'(A80)')  parfilein(i)  ! input CTF parameter information
  READ (5,'(A80)')  parfileout(i)  ! output parameter file
END DO

! DATA STREAMS:
! 1 = Particles
! 2 = Map
! 3 = Parameter file input
! 5 = RESERVED FOR EXECUTING SCRIPT
! 4 = Parameter file output
! 6 = RESERVED FOR SCREEN OUTPUT
! 7 = REMOVED: Orientation list file (output)

!-----------------------------------------
! Open map 
!-----------------------------------------
CALL Imopen(2,inmap,"RO")
CALL Irdhdr(2,nxyz,mxyz,mode,dmin,dmax,dmean)
nx         = nxyz(1)
ny         = nxyz(2)
nz         = nxyz(3)
nparts = (ilast-ifirst)+1
nprojs = (ncheck/2 * INT(90/psi_step))*8
CALL GenOptimizationConstants

! Restate resolution limits in terms of pixel radii
rmax1=(psize*nx)/rmax1
IF (rmax1<1) rmax1=1.01 ! Never use the origin of the FT for alignment
rmax1sq=rmax1**2
rmax2=(psize*nx)/rmax2
rmax2sq=rmax2**2
ri=ri/psize
nparts = (ilast-ifirst)+1
psi_step=psi_step*pi/180

!-----------------------------------------
! Format statements
!-----------------------------------------
format1='(I7,5F8.3,F8.0,I6,2F9.1,F8.2,F7.2)' ! Parameter file format
format6='(I7,5F8.3,F8.0,I6,2F9.1,F8.2,F7.2,I3.2,A1,I2.2,A1,I2.2)' ! Parameter with time stamp
format2='(A20,I7,A6,I3.2,A1,I2.2,A1,I2.2)'
format3='(I7,2F9.3)'   
format4='(F8.2,2I3,F6.3)'
format5='(A9,I6)'

!-------------------------------------------------
! Allocate arrays
!-------------------------------------------------
radius = 0.5 * nx
! Volumes images
ALLOCATE (map(nx,ny,nz))                     ; map=0
ALLOCATE (workingvol(nxby2plus1,ny,nz))      ; workingvol=0

! Images
ALLOCATE (image(nx,ny))                  ; image=0
ALLOCATE (imagec(nxby2plus1,ny))         ; imagec=0
ALLOCATE (workingslice(nxby2plus1,ny))   ; workingslice=0
ALLOCATE (workingslice2(nxby2plus1,ny))  ; workingslice2=0
ALLOCATE (workingslice3(nxby2plus1,ny))  ; workingslice3=0
ALLOCATE (ctfx(nxby2plus1,ny))           ; ctfx=0
ALLOCATE (ctfx_wrap(nxby2plus1,ny))      ; ctfx_wrap=0
ALLOCATE (particle(nx,ny))               ; particle=0
ALLOCATE (particlec(nxby2plus1,ny))      ; particlec=0
ALLOCATE (ccfmask(nx,ny))                ; ccfmask=0
ALLOCATE (rimask(nx,ny))                 ; rimask=0

! Normalization
ALLOCATE (sumsquaredref(nprojs))         ; sumsquaredref=0

! FFTW variables
ALLOCATE (image_fftw(nx,ny))             ; image_fftw=0
ALLOCATE (imagec_fftw(nxby2plus1,ny))    ; imagec_fftw=0
ALLOCATE (vol_fftw(nx,ny,nz))            ; vol_fftw=0
ALLOCATE (volc_fftw(nxby2plus1,ny,nz))   ; volc_fftw=0

!-------------------------------------------------
! Plan Fourier Transforms
!-------------------------------------------------
CALL Dfftw_plan_dft_r2c_2d(plan2df,nx,ny,image_fftw,imagec_fftw,fftw_measure)
CALL Dfftw_plan_dft_c2r_2d(plan2dr,nx,ny,imagec_fftw,image_fftw,fftw_measure)
CALL Dfftw_plan_dft_r2c_3d(plan3df,nx,ny,nz,vol_fftw,volc_fftw,fftw_estimate)
CALL Dfftw_plan_dft_c2r_3d(plan3dr,nx,ny,nz,volc_fftw,vol_fftw,fftw_estimate)
!-----------------------------------------
! Prepare CTF constants
!-----------------------------------------
AKV=1000.0*AKV
WL=12.3/SQRT(AKV+AKV**2/(10.0**6.0))
wgh1=SQRT(1.0-wgh**2)
wgh2=wgh
thetatr=wl/(psize*nx)

!-------------------------------------------------
! Prepare a mask for the CCF
!-------------------------------------------------
WRITE (*,*) "Preparing CCF and RI masks"
DO lin=1,ny
  DO vox=1,nx
    IF (lin<=nyby2plus1) ypos = lin-1      
    IF (lin>nyby2plus1)  ypos = lin-ny-1
    IF (vox<=nxby2)   xpos = vox-1      
    IF (vox>nxby2)    xpos = vox-nx-1
    radius=(xpos**2 + ypos**2)**0.5
    IF (radius<=shiftmax) THEN
      ccfmask(vox,lin)=.true.
    ELSE
      ccfmask(vox,lin)=.false.
    END IF
  END DO
END DO


!-------------------------------------------------
! Prepare a mask for the real space image
!-------------------------------------------------
DO lin=1,ny
  DO vox=1,nx
    xpos=vox-nxby2
    ypos=lin-nyby2
    radius=(xpos**2 + ypos**2)**0.5
    IF (radius<=ri) THEN
      rimask(vox,lin)=1
    ELSE
      rimask(vox,lin)=0
    END IF
  END DO
END DO

!-------------------------------------------------
! Read in and Fourier transform map
!-------------------------------------------------
DO sec=1,nz
  DO lin=1,ny
    CALL Imposn(2,sec-1,lin-1) 
    CALL Irdlin (2,map(:,lin,sec),*999)
  END DO 
END DO 
vol_fftw=DBLE(map)
CALL Dfftw_execute(plan3df)

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
! Get to correct starting point for stack, parameter file, and CTF file
!-------------------------------------------------
! Open first particle stack
CALL Imopen(1,instack(1),"RO")
CALL Irdhdr(1,nxyz,mxyz,mode,dmin,dmax,dmean)
IF (nx.ne.nxyz(1)) STOP 'particle stack and map do not have matching dimensions'
IF (ny.ne.nxyz(2)) STOP 'particle stack and map do not have matching dimensions'
stackparts = nxyz(3)
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
OPEN(unit=3,file=parfilein(stack),status="old")
! Skip comment lines
DO
  READ (3,'(A100)') line
  IF (line(1:23).eq."C           PSI   THETA") EXIT
END DO

!-----------------------------------------
! Open parameter file output
!-----------------------------------------
OPEN(unit=4,file=parfileout(stack),status="unknown")
WRITE (4,'(A30)')  "C  Search-Fspace parameter file" 
WRITE (4,'(A1)')   "C"
WRITE (4,'(A94)')  "C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMAX"


!-------------------------------------------------
! Main Loop over all particles
!-------------------------------------------------
WRITE (*,*) "Main loop"
stackposition=stackposition-1
Eachpart: DO part=ifirst,ilast
  stackposition=stackposition+1
  CALL Itime(now)
  !-------------------------------------------------
  ! Open a new stack and parameter file if necessary
  !-------------------------------------------------
  IF (part>described) THEN
    write (*,*) "opening new stack and parameter file",described
    stack=stack+1
    CALL Imclose(1)
    CALL Imopen(1,instack(stack),"RO")
    CALL Irdhdr(1,nxyz,mxyz,mode,dmin,dmax,dmean)
    stackparts = nxyz(3)
    stackposition=1
    described=described+stackparts
    CLOSE(3) ! Input parfile
    CLOSE(4) ! Output parfile
    OPEN(unit=3,file=parfilein(stack),status="unknown")
    OPEN(unit=4,file=parfileout(stack),status="unknown")
    DO
      READ (3,'(A100)') line
      IF (line(1:23).eq."C           PSI   THETA") EXIT
    END DO
    ! Prepare new output parfile
    WRITE (4,'(A30)')  "C  Align-Fspace parameter file" 
    WRITE (4,'(A1)')   "C"
    WRITE (4,'(A94)')  "C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMAX"
  END IF
  !-------------------------------------------------
  ! read in the relevant line of the parameter file
  !-------------------------------------------------
  DO
    READ(3,format1) ilist,psibest,thetabest,phibest,shxbest,shybest,&
                    amagp,film,dfmid1,dfmid2,angasti,jnkr
    IF (ilist==stackposition) EXIT
  END DO
  angast=angasti*piby180
  psibest=0
  thetabest=0
  phibest=0
  shxbest=0
  shybest=0
  CCbest=0
  !!shxbest=2*pi*shx/nx
  !!shybest=2*pi*shy/ny
  !-------------------------------------------------
  ! Generate a new CTF multiplier if necessary (not wraparound order)
  !-------------------------------------------------
  IF (dfmid1last.ne.dfmid1.or.dfmid2last.ne.dfmid2.or.angastlast.ne.angasti) THEN
    WRITE (*,'(A19,I8)') "New CTF at particle",part
    dfmid1last=dfmid1
    dfmid2last=dfmid2
    angastlast=angasti
    ctfx=0
    CALL GenCTFmultiplier(ctfx,ctfx_wrap,nxby2,nyby2,ny,rmax1sq,rmax2sq,cs,&
          wl,wgh1,wgh2,dfmid1,dfmid2,angast,thetatr)
    newctf=.true.
  ELSE
    newctf=.false.
  END IF
  !-------------------------------------------------
  ! put image values into array image(:,:)
  !-------------------------------------------------
  CALL Imposn(1,stackposition-1,0)
  CALL Irdsec (1,image(:,:),*999)    
  particle=DBLE(image)
  !-------------------------------------------------
  ! Transform particle image and apply CTF
  !-------------------------------------------------
  ! NB: for speedup, particlec is now the CONJG of the particle transform
  image_fftw=DBLE(particle)
  CALL Dfftw_execute(plan2df)
  particlec=CONJG(imagec_fftw*ctfx_wrap**ctfexppart)
  !-------------------------------------------------
  ! Get sumsquaredpart for normalized CCF calculation later and apply mask
  !-------------------------------------------------
  imagec_fftw(:,:)=particlec(:,:)
  CALL Dfftw_execute(plan2dr)
  image_fftw=image_fftw*rimask*nxnyinv
  particle(:,:)=image_fftw
  sumsquaredpart=0
  sumsquaredpart=SUM(REAL(particle)**2)
  CALL Dfftw_execute(plan2df) 
  particlec=imagec_fftw
  !-------------------------------------------------
  ! Cycle through all Euler angles
  !-------------------------------------------------
  orientationID=0
  CCBest=0
  Eachncheck: DO i=1,ncheck/2
    h=-1+2*((i-1)/(REAL(ncheck)-1))
    theta=acos(h)
    IF(i==1.or.i==ncheck) THEN 
      phi=0
    ELSE
      phi= amod((lastphi+3.6/sqrt(ncheck*(1-h**2))),2*pi)
    END IF
    lastphi=phi
    !-------------------------------------------------
    ! Search all psi
    !-------------------------------------------------
    Eachpsi: DO psi=0,pi/2-(1E-5),psi_step
      orientationID=orientationID+1
      !-------------------------------------------------
      ! Extract Fourier slice at given orientation
      !-------------------------------------------------
      CALL Extractimage_noshifts(workingvol,workingslice,nx,ny,nz,phi,theta,psi,rmax1sq,rmax2sq)
      !-------------------------------------------------
      ! Apply the CTF multiplier to the extracted slice
      !-------------------------------------------------
      workingslice=workingslice*CMPLX(ctfx)**ctfexpmod
      !-------------------------------------------------
      ! For a given extracted slice, quickly check the opposite hemisphere and 90' rotations
      ! keep track of these transformations on psi, theta and phi
      !-------------------------------------------------
      ! Check rotations of the projection by pi/2, pi and 3pi/2
      Eachquadrant: DO quadrant=1,4
        IF (quadrant==1) THEN
          psiout=psi
          workingslice2=workingslice
        ELSE IF (quadrant==2) THEN
          psiout=psi+pi/2
          IF (psiout>2*pi) psiout=psiout-2*pi
          DO k=1,nyby2
            DO j=1,nxby2plus1
               workingslice2(j,k)=workingslice(nyby2plus1-k,nxby2-1+j)
            END DO
          END DO
          DO k=nyby2plus1,ny
            DO j=1,nxby2
              workingslice2(j,k)=CONJG(workingslice(k+1-nyby2,nxby2plus1-j))
            END DO
            workingslice2(nxby2plus1,k)=CONJG(workingslice(k+1-nyby2,nx+1))
          END DO
        ELSE IF (quadrant==3) THEN
          psiout=psi+pi
          IF (psiout>2*pi) psiout=psiout-2*pi
          workingslice2=CONJG(workingslice)
        ELSE IF (quadrant==4) THEN
          psiout=psi+3*pi/2
          IF (psiout>2*pi) psiout=psiout-2*pi
          DO k=1,nyby2-1
            DO j=1,nxby2plus1
              workingslice2(j,k)=CONJG(workingslice(nyby2plus1-k,nxby2-1+j))
            END DO
          END DO
          DO k=nyby2,ny
            DO j=1,nxby2
              workingslice2(j,k)=workingslice(k+1-nyby2,nxby2plus1-j)
            END DO
            workingslice2(nxby2plus1,k)=workingslice(k+1-nyby2,nx+1)
          END DO
        END IF
        ! Check mirrors of the projection
        Eachhemisphere: DO hemisphere=1,2
          IF (hemisphere==1) THEN
           ! Case 1: The projection is used as is
            phiout=phi
            thetaout=theta
            psiout=psiout
            workingslice3=workingslice2
          ELSE
            ! Case 2: The mirror image of the projection is used
            ! adjust Euler angles to reflect mirroring
            ! phi=phi+pi; theta=pi-theta; psi=-psi
            workingslice3(:,ny)=workingslice2(:,ny)
            DO j=1,ny-1
              workingslice3(:,j)=workingslice2(:,ny-j)
            END DO
            phiout=phi+pi
            IF (phiout.gt.2*pi) phiout=phiout-2*pi
            thetaout=pi-theta
            psiout=-psiout
            IF (psiout.lt.0) psiout=psiout+2*pi
          END IF
          ! Unload working slice into wraparound order
          imagec(1:nxby2plus1,nyby2plus2:ny)=workingslice3(1:nxby2plus1,1:nyby2minus1)
          imagec(1:nxby2plus1,1:nyby2plus1)=workingslice3(1:nxby2plus1,nyby2:ny)
          ! Move origin
          CALL Moveorigin2d(imagec,nxby2,nxby2plus1,ny)
          !-----------------------------------------------------------------
          ! If this is a new CTF AND hemisphere=1 and quadrant=1 get a sum of squared pixel value 
          ! It is unclear to me why this sum-squared value varies so much for each psi (but it does)
          !-----------------------------------------------------------------
          IF (newctf.and.hemisphere==1.and.quadrant==1)  THEN
            !write (*,*) "recalculating sumsquaredref"
            ! FT Into realspace
            imagec_fftw=imagec
            CALL Dfftw_execute(plan2dr)
            image_fftw=image_fftw*nxnyinv
            ! Calculate sum-squared pixel value of reference
            ssref=0
            ssref=SUM(REAL(image_fftw)**2)
            sumsquaredref(orientationID)=ssref
          END IF 
          !imagec=workingslice3
          !-----------------------------------------------------------------
          ! Calculate normalized CCF
          !-----------------------------------------------------------------
          ! NB should be able to replace above with
          imagec_fftw=imagec*particlec
          CALL Dfftw_execute(plan2dr)
          ccnorm=1/SQRT(sumsquaredref(orientationID)*sumsquaredpart)
          CCC=MAXVAL(image_fftw,ccfmask)*nxnyinv*ccnorm  ! Normalize and only look within mask
          IF (CCC>CCbest) THEN
            CCBest=CCC
            psibest=psiout
            thetabest=thetaout
            phibest=phiout
            shifts=MAXLOC(image_fftw,ccfmask)
            shx=shifts(1)
            shy=shifts(2)
            IF (shy<nyby2plus1)  shybest=REAL(-1*shy+1)
            IF (shy>nyby2plus1) shybest=REAL(ny-shy+1)
            IF (shx<nxby2plus1) shxbest=REAL(-1*shx+1)
            IF (shx>nxby2plus1)  shxbest=REAL(nx-shx+1)
          END IF
        END DO Eachhemisphere
      END DO Eachquadrant
    END DO Eachpsi
    !-----------------------------------------------------------------
    ! Output CC, psi, SHX, SHY at each ncheck
    !-----------------------------------------------------------------
  END DO Eachncheck
  ! Output best choice here
  WRITE (*,format6) part,psibest*oneeightybypi,thetabest*oneeightybypi,phibest*oneeightybypi,shxbest,&
                 shybest,amagp,film,dfmid1,dfmid2,angast*180/pi,CCbest*100,now(1),":",now(2),":",now(3)
  WRITE (4,format1) part,psibest*oneeightybypi,thetabest*oneeightybypi,phibest*oneeightybypi,shxbest,&
                 shybest,amagp,film,dfmid1,dfmid2,angast*180/pi,CCbest*100
END DO Eachpart

!-------------------------------------------------
! Close files and destroy plans
!-------------------------------------------------
CALL Imclose(1)
CALL Imclose(2)
CLOSE(3)
CLOSE(4)

CALL Dfftw_destroy_plan(plan2df)
CALL Dfftw_destroy_plan(plan2dr)
CALL Dfftw_destroy_plan(plan3df)
CALL Dfftw_destroy_plan(plan3dr)

STOP "Normal termination of Search_fspace"
999 STOP "End of File Read Error"
END PROGRAM Search_fspace
