! JLR ! Version 1.01 Dated 171111 
! 10/12/10 First separatin of subroutines and functions from programs
! 22/04/11 Conversion to modular modern Fortran

MODULE Consts
  REAL, PARAMETER :: pi   =3.1415926535897
  REAL, PARAMETER :: twopi=6.2831853071796
  REAL, PARAMETER :: piby180=pi/180
END MODULE Consts

MODULE Imageinfo
  INTEGER          :: nx,ny,nz
  INTEGER          :: nxby2,nyby2,nzby2
  INTEGER          :: nxby2plus1,nyby2plus1,nzby2plus1,nyby2minus1
  INTEGER          :: nyby2plus2
  DOUBLE PRECISION :: nxnyinv,nxnynzinv
CONTAINS
  SUBROUTINE GenOptimizationConstants
  IMPLICIT NONE
  nxby2=nx/2            ; nyby2=ny/2         ; nzby2=nz/2
  nxby2plus1=nxby2+1    ; nyby2plus1=nyby2+1 ; nzby2plus1=nzby2 + 1
  nyby2minus1=nyby2-1   ; nyby2plus2 = nyby2+2
  nxnyinv=1/DBLE(nx*ny) ; nxnynzinv=1/DBLE(nx*ny*nz)
  END SUBROUTINE GenOptimizationConstants
END MODULE

MODULE CTFmod
CONTAINS
  !-------------------------------------------------------
  SUBROUTINE GenCTFmultiplier(ctfx,ctfx_wrap,nxby2,nyby2,ny,rmax1sq,rmax2sq,CS,WL,WGH1,WGH2,DFMID1,DFMID2,ANGAST,THETATR)
  !------------------------------------------------------
  IMPLICIT NONE
  INTEGER          :: kxarray,kyarray,vox,lin,nxby2,nyby2,ny,kxwrap,kywrap
  REAL             :: radiussq,rmax1sq,rmax2sq
  REAL             :: cs,wl,wgh1,wgh2,dfmid1,dfmid2,angast,thetatr!,ctf
  DOUBLE PRECISION :: ctfx(nxby2+1,ny), ctfx_wrap(nxby2+1,ny) 
  ctfx=0
  ctfx_wrap=0
  DO kyarray=-nyby2+1,nyby2,1 
    DO kxarray=0,nxby2,1
      vox = kxarray+1
      lin = kyarray+nyby2                  
      IF (lin<nyby2)  kywrap=lin+nyby2+1
      IF (lin>=nyby2) kywrap=lin-nyby2+1
      radiussq=kxarray**2+kyarray**2
      IF(radiussq>=rmax1sq.and.radiussq<=rmax2sq) THEN 
        ctfx(vox,lin)=DBLE(CTF(CS,WL,WGH1,WGH2,DFMID1,DFMID2,ANGAST,THETATR,kxarray,kyarray))
        ctfx_wrap(vox,kywrap)=ctfx(vox,lin)
      ELSE
        ctfx(vox,lin)=0
        ctfx_wrap(vox,kywrap)=0
      END IF
    END DO
  END DO
  RETURN
  END SUBROUTINE GenCTFMultiplier
  !-------------------------------------------------------
  REAL FUNCTION CTF(CSf,WLf,WGH1f,WGH2f,dfmid1f,dfmid2f,angastf,THETATRf,IX,IY)
  !-------------------------------------------------------
  !  Calculates contrast transfer function, including the contribution from 
  !     amplitude contrast WGH2 - WGH1 is the resulting phase contrast
  !
  USE Consts
  IMPLICIT NONE
  REAL       csf,wlf,wgh1f,wgh2f,dfmid1f,dfmid2f,angastf,thetatrf 
  INTEGER    ix,iy
  REAL       rad,angle,angspt,c1,c2,angdif,ccos,df,chi
  rad=ix**2+iy**2
  IF (rad.NE.0.0) THEN
    rad=SQRT(rad)
    angle=rad*THETATRf
    angspt=ATAN2(REAL(iy),REAL(ix))
    c1=twopi*angle*angle/(2.0*wlf)
    c2=-c1*csf*angle*angle/2.0
    angdif=angspt-angastf
    ccos=COS(2.0*angdif)
    DF=0.5*(dfmid1f+dfmid2f+ccos*(dfmid1f-dfmid2f))
    chi=c1*df+c2
    ctf=-wgh1f*SIN(chi)-wgh2f*COS(chi)
  ELSE
    ctf=-wgh2f
  ENDIF
  END FUNCTION CTF
END MODULE CTFmod

MODULE Fouriervolume
  USE Consts
  IMPLICIT NONE
  REAL, PARAMETER, PRIVATE          :: irad = 1 !! Determines amount of interpolation
  REAL, PRIVATE                     :: kx,ky,kz,kxplane,kyplane,kzplane,radius,radiussq
  INTEGER, PRIVATE                  :: kxs,kxf,kys,kyf,kzs,kzf,kxvol,kyvol,kzvol
  INTEGER, PRIVATE                  :: kxarray,kyarray,kzarray,kxwrap,kywrap,kzwrap,kxarraysq
  REAL, DIMENSION(6), PRIVATE       :: extract,insert
  REAL, PRIVATE                     :: cpsi,spsi,cthe,sthe,cphi,sphi
  REAL, PRIVATE                     :: boxftv,phase,sampctf,shftx,shfty
  INTEGER, PRIVATE                  :: hermitian
  COMPLEX, PRIVATE                  :: samp, pshft
  INTEGER, PRIVATE                  :: nxby2,nyby2,nzby2,i,j
  ! The following line generates the sinclut. Initialized only once during execution.
  REAL, PARAMETER, PRIVATE, DIMENSION(2000) :: &
    sinclut=[(SIN(FLOAT(j)*piby180)/(FLOAT(j)*piby180),j=1,2000)]
CONTAINS
  !-------------------------------------------------------
  SUBROUTINE Addimage(workingvol,normvol,workingslice,nx,ny,nz,phi,theta,psi,shx,shy,rmax2sq,ctfx,ctfexppart)
  !-------------------------------------------------------
  ! Insert an image into the working volume
  ! Receives workingvol, normvol, workingslice,ctfx (all not in wraparound order)
  !          nx,ny,nz (all in pixels)
  !          phi,theta,psi (all in radians)
  !          shx,shy (in pixels)
  !          rmax2sq: pixel radius to include in F-space squared
  ! JLR - partially optimized
  USE Consts
  IMPLICIT NONE
  ! Calling variables
  REAL, INTENT(IN)                :: psi,theta,phi,shx,shy
  INTEGER, INTENT(IN)             :: nx,ny,nz
  DOUBLE PRECISION, INTENT(IN)    :: ctfx(nx/2+1,ny)
  DOUBLE PRECISION, INTENT(INOUT) :: normvol(nx/2+1,ny,nz)
  DOUBLE COMPLEX, INTENT (INOUT)  :: workingvol(nx/2+1,ny,nz)
  DOUBLE COMPLEX, INTENT (IN)     :: workingslice(nx/2+1,ny)
  REAL, INTENT(IN)                :: ctfexppart, rmax2sq
  ! Other variable
  REAL                            :: arg(3)

  shftx=twopi*shx/nx ! Multiply shifts by 2pi/nx
  shfty=twopi*shy/ny

  nxby2=nx/2
  nyby2=ny/2
  nzby2=nz/2

  cphi=COS(phi)
  sphi=SIN(phi)
  cthe=COS(theta)
  sthe=SIN(theta)
  cpsi=COS(psi)
  spsi=SIN(psi)
  insert(1)=cphi*cthe*cpsi-sphi*spsi
  insert(2)=sphi*cthe*cpsi+cphi*spsi
  insert(3)=-sthe*cpsi
  insert(4)=-cphi*cthe*spsi-sphi*cpsi
  insert(5)=-sphi*cthe*spsi+cphi*cpsi
  insert(6)=sthe*spsi
  ! determine where pixels in the working slice map onto the working volume after application of reversed Euler angles
  Eachkx: DO kxarray=0,nxby2,1
    Eachky: DO kyarray=-nyby2+1,nyby2,1 ! Do not need to insert ky=-fc because that is the same as ky=fc
      ! Determine value from 2D array to be deposited in 3D array
      phase=shftx*kxarray+shfty*kyarray
      pshft=CMPLX(COS(phase),SIN(phase))
      samp=CMPLX(workingslice(kxarray+1,kyarray+ny/2))*pshft
      sampctf=(ctfx(kxarray+1,kyarray+nyby2))**ctfexppart ! Here CTFx needs to be out of wraparound order
      radiussq=kxarray**2+kyarray**2
      IF(radiussq<=rmax2sq) THEN ! only work out to rmax2      
        ! Determine location in 3D array where value should be placed (off lattice)
        kxplane = insert(1)*(kxarray)+insert(4)*(kyarray) 
        kyplane = insert(2)*(kxarray)+insert(5)*(kyarray)
        kzplane = insert(3)*(kxarray)+insert(6)*(kyarray)
        ! shift from fractional working volume coordinates to "fictional"fractional array coordinates
        kxplane = kxplane+nxby2+1   
        kyplane = kyplane+nyby2+1
        kzplane = kzplane+nzby2+1
        ! Interpolation
        ! assign ranges for voxels to be considered
        kxs = INT(kxplane)-irad+1
        kxf = INT(kxplane)+irad
        kys = INT(kyplane)-irad+1
        kyf = INT(kyplane)+irad
        kzs = INT(kzplane)-irad+1
        kzf = INT(kzplane)+irad
        ! Scan over voxels in workingvol that will receive values from off-lattice pixels in image
        Intkx: DO kxvol=kxs,kxf
          Intky: DO kyvol=kys,kyf
            Intkz: DO kzvol=kzs,kzf
              kxwrap=kxvol-nxby2 ! Actual array positions (x can be negative, y and z can be 0)
              kywrap=kyvol-1
              kzwrap=kzvol-1
              IF (kxwrap<1) THEN     ! If a point is found in the hemisphere not covered (but Hermitian)
                kxwrap=-1*kxwrap+2
                kywrap=ny-kywrap               
                kzwrap=nz-kzwrap
                hermitian=-1
              ELSE
                hermitian=1
              END IF
              IF (kywrap==0) kywrap=ny ! take care of Nyquist voxels
              IF (kzwrap==0) kzwrap=nz
              IF (kxwrap>=1.and.kxwrap<=nxby2+1.and.kywrap>=1.and.kywrap<=ny.and.kzwrap>=1.and.kzwrap<=nz) THEN
                !IF (kxvol>=1.and.kxvol<=nx+1.and.kyvol>=0.and.kyvol<=ny.and.kzvol>=0.and.kzvol<=nz) THEN
                arg(1)=(kxplane-kxvol)
                arg(2)=(kyplane-kyvol)
                arg(3)=(kzplane-kzvol)
                boxftv=boxft_lut(arg)
                IF (hermitian==1) THEN
                  workingvol(kxwrap,kywrap,kzwrap)=workingvol(kxwrap,kywrap,kzwrap)+samp*sampctf*boxftv*boxftv
                ELSE IF (hermitian==-1) THEN
                  workingvol(kxwrap,kywrap,kzwrap)=workingvol(kxwrap,kywrap,kzwrap)+CONJG(samp)*sampctf*boxftv*boxftv
                END IF
                normvol(kxwrap,kywrap,kzwrap)   =normvol(kxwrap,kywrap,kzwrap)   +(sampctf*boxftv)**2
              END IF
            END DO Intkz
          END DO Intky
        END DO Intkx
      END IF
    END DO Eachky
  END DO Eachkx
  RETURN
  END SUBROUTINE Addimage

  !-------------------------------------------------------
  SUBROUTINE Extractimage(workingvol,workingslice,nx,ny,nz,phi,theta,psi,shx,shy,rmax1sq,rmax2sq)
  !-------------------------------------------------------
  ! Extract Fourier slice at given orientation, angles input in radians, shifts in pixels
  ! JLR (some optimizations from CM)
  ! JLR Shifts allowed: input as pixel shifts, corrected within subroutine
  USE Consts
  IMPLICIT NONE
  ! Calling variables
  REAL, INTENT(IN)            :: psi,theta,phi,shx,shy
  INTEGER, INTENT(IN)         :: nx,ny,nz
  DOUBLE COMPLEX, INTENT(OUT) :: workingvol(nx/2+1,ny,nz)
  DOUBLE COMPLEX, INTENT(OUT) :: workingslice(nx/2+1,ny)
  REAL, INTENT(IN)            :: rmax1sq,rmax2sq
  ! Other
  REAL                        :: arg(3)

  nxby2=nx/2
  nyby2=ny/2
  nzby2=nz/2
  shftx=-shx*2*pi/nx
  shfty=-shy*2*pi/ny
  workingslice=0
  ! Multiply the combined Eulerian rotation matrix (with the reversed rotations) by the vector [x,y,0]
  cphi=COS(phi)
  sphi=SIN(phi)
  cthe=COS(theta)
  sthe=SIN(theta)
  cpsi=COS(psi)
  spsi=SIN(psi)
  extract(1)=cphi*cthe*cpsi-sphi*spsi
  extract(2)=sphi*cthe*cpsi+cphi*spsi
  extract(3)=-sthe*cpsi
  extract(4)=-cphi*cthe*spsi-sphi*cpsi
  extract(5)=-sphi*cthe*spsi+cphi*cpsi
  extract(6)=sthe*spsi
  ! determine where pixels in the working slice map onto the working volume after application of reversed Euler angles
  Eachkx: DO kxarray=0,nxby2,1
    kxarraysq=kxarray*kxarray
    IF (kxarraysq>rmax2sq) EXIT
    Eachky: DO kyarray=-nyby2+1,nyby2,1 ! Do not need to obtain ky=-fc because that is the same as ky=fc
      radiussq=(kxarraysq+kyarray*kyarray)  ! NB: we are now working with radiussq, kxarraysq, and kyarraysq
      IF(.not.(radiussq>rmax2sq.or.radiussq<rmax1sq)) THEN ! only work within the ring between rmax1 and rmax2
        phase=shftx*kxarray+shfty*kyarray  
        pshft=CMPLX(COS(phase),SIN(phase))
        !pshft=CMPLX(1,0) !!!!!!!!!!!!1   Forces shifts to 0
        kxplane = extract(1)*(kxarray)+extract(4)*(kyarray) 
        kyplane = extract(2)*(kxarray)+extract(5)*(kyarray)
        kzplane = extract(3)*(kxarray)+extract(6)*(kyarray)
        ! shift from fractional working volume coordinates to fractional array coordinates
        kxplane = kxplane+nxby2+1  ! "fictional" always positive array positions
        kyplane = kyplane+nyby2+1
        kzplane = kzplane+nzby2+1
        ! Interpolation
        ! assign ranges for voxels to be considered
        kxs = INT(kxplane)-irad+1
        kxf = INT(kxplane)+irad
        kys = INT(kyplane)-irad+1
        kyf = INT(kyplane)+irad
        kzs = INT(kzplane)-irad+1
        kzf = INT(kzplane)+irad
        samp=CMPLX(0,0)
        ! Scan over voxels in workingvol that will contribute to off-lattice pixels in workingslice
        Intkx: DO kxvol=kxs,kxf
          Intky: DO kyvol=kys,kyf
            Intkz: DO kzvol=kzs,kzf
              kxwrap=kxvol-nxby2 ! Actual array positions (x can be negative, y and z can be 0)
              kywrap=kyvol-1
              kzwrap=kzvol-1
              IF (kxwrap<1) THEN     ! If a point is found in the hemisphere not covered (but Hermitian)
                kxwrap=-1*kxwrap+2
                kywrap=ny-kywrap               
                kzwrap=nz-kzwrap
                hermitian=-1
              ELSE
                hermitian=1
              END IF
              IF (kywrap==0) kywrap=ny
              IF (kzwrap==0) kzwrap=nz
              IF (kxwrap>=1.and.kxwrap<=nxby2+1.and.kywrap>=1.and.kywrap<=ny.and.kzwrap>=1.and.kzwrap<=nz) THEN
                arg(1)=(kxplane-kxvol)
                arg(2)=(kyplane-kyvol)
                arg(3)=(kzplane-kzvol)
                boxftv=boxft_lut(arg)
                IF (hermitian==1) THEN
                  samp=samp+workingvol(kxwrap,kywrap,kzwrap)*boxftv
                ELSE IF (hermitian==-1) THEN
                  samp=samp+CONJG(workingvol(kxwrap,kywrap,kzwrap))*boxftv
                END IF
              ELSE
                samp=samp
              END IF
            END DO Intkz
          END DO Intky
        END DO Intkx
        ! Take the interpolated value and put into working slice, apply phase shift
        workingslice(kxarray+1,kyarray+nyby2)=CMPLX(samp)*pshft
      END IF
    END DO Eachky
  END DO Eachkx
  RETURN
  END SUBROUTINE Extractimage

  !-------------------------------------------------------
  SUBROUTINE Extractimage_noshifts(workingvol,workingslice,nx,ny,nz,phi,theta,psi,rmax1sq,rmax2sq)
  !-------------------------------------------------------
  ! Extract Fourier slice at given orientation
  ! JLR (some optimizations from CM)
  ! JLR Removed phase shift to account for usage in this context
  IMPLICIT NONE
  ! Calling variables
  REAL, INTENT(IN)               :: psi,theta,phi
  INTEGER, INTENT(IN)            :: nx,ny,nz
  DOUBLE COMPLEX, INTENT(IN)     :: workingvol(nx/2+1,ny,nz)
  DOUBLE COMPLEX, INTENT(OUT)    :: workingslice(nx/2+1,ny)
  REAL, INTENT(IN)               :: rmax1sq,rmax2sq
  ! Other variable
  REAL                            :: arg(3)

  nxby2=nx/2
  nyby2=ny/2
  nzby2=nz/2
  workingslice=0

  ! Multiply the combined Eulerian rotation matrix (with the reversed rotations) by the vector [x,y,0]
  cphi=COS(phi)
  sphi=SIN(phi)
  cthe=COS(theta)
  sthe=SIN(theta)
  cpsi=COS(psi)
  spsi=SIN(psi)
  extract(1)=cphi*cthe*cpsi-sphi*spsi
  extract(2)=sphi*cthe*cpsi+cphi*spsi
  extract(3)=-sthe*cpsi
  extract(4)=-cphi*cthe*spsi-sphi*cpsi
  extract(5)=-sphi*cthe*spsi+cphi*cpsi
  extract(6)=sthe*spsi

  ! determine where pixels in the working slice map onto the working volume after application of reversed Euler angles
  Eachkx: DO kxarray=0,nxby2,1
    kxarraysq=kxarray*kxarray
    IF (kxarraysq>rmax2sq) EXIT
    Eachky: DO kyarray=-nyby2+1,nyby2,1 ! Do not need to obtain ky=-fc because that is the same as ky=fc
      radiussq=(kxarraysq+kyarray*kyarray)  ! NB: we are now working with radiussq, kxarraysq, and kyarraysq
      IF(.not.(radiussq>rmax2sq.or.radiussq<rmax1sq)) THEN ! only work within the ring between rmax1 and rmax2
        pshft=CMPLX(1,0)
        kxplane = extract(1)*(kxarray)+extract(4)*(kyarray) 
        kyplane = extract(2)*(kxarray)+extract(5)*(kyarray)
        kzplane = extract(3)*(kxarray)+extract(6)*(kyarray)
        ! shift from fractional working volume coordinates to fractional array coordinates
        kxplane = kxplane+nxby2+1  ! "fictional" always positive array positions
        kyplane = kyplane+nyby2+1
        kzplane = kzplane+nzby2+1
        ! Interpolation
        ! assign ranges for voxels to be considered
        kxs = INT(kxplane)-irad+1
        kxf = INT(kxplane)+irad
        kys = INT(kyplane)-irad+1
        kyf = INT(kyplane)+irad
        kzs = INT(kzplane)-irad+1
        kzf = INT(kzplane)+irad
        samp=CMPLX(0,0)
        ! Scan over voxels in workingvol that will contribute to off-lattice pixels in workingslice
        Intkx: DO kxvol=kxs,kxf
          Intky: DO kyvol=kys,kyf
            Intkz: DO kzvol=kzs,kzf
              kxwrap=kxvol-nxby2 ! Actual array positions (x can be negative, y and z can be 0)
              kywrap=kyvol-1
              kzwrap=kzvol-1
              IF (kxwrap<1) THEN     ! If a point is found in the hemisphere not covered (but Hermitian)
                kxwrap=-1*kxwrap+2
                kywrap=ny-kywrap               
                kzwrap=nz-kzwrap
                hermitian=-1
              ELSE
                hermitian=1
              END IF
              IF (kywrap==0) kywrap=ny
              IF (kzwrap==0) kzwrap=nz
              IF (kxwrap>=1.and.kxwrap<=nxby2+1.and.kywrap>=1.and.kywrap<=ny.and.kzwrap>=1.and.kzwrap<=nz) THEN
                arg(1)=(kxplane-kxvol)
                arg(2)=(kyplane-kyvol)
                arg(3)=(kzplane-kzvol)
                boxftv=boxft_lut(arg)
                IF (hermitian==1) THEN
                  samp=samp+workingvol(kxwrap,kywrap,kzwrap)*boxftv
                ELSE IF (hermitian==-1) THEN
                  samp=samp+CONJG(workingvol(kxwrap,kywrap,kzwrap))*boxftv
                END IF
              ELSE
                samp=samp
              END IF
            END DO Intkz
          END DO Intky
        END DO Intkx
        ! Take the interpolated value and put into working slice, apply phase shift
        workingslice(kxarray+1,kyarray+nyby2)=CMPLX(samp)*pshft
      END IF
    END DO Eachky
  END DO Eachkx
  RETURN
  END SUBROUTINE Extractimage_noshifts

  !-------------------------------------------------------
  SUBROUTINE Renormalizemap(workingvol,normingvol,nx,ny,nz)
  !-------------------------------------------------------
  ! JLR: Renormalize the working volume
  !      9/04/09  Modified to overwrite workingvolume to save memory usage
  !               Norming volume converted to DP from Double Complex
  !      10/08/10 Fixed bug: set norm to 0 before calcuating it
  ! Maps not in wraparound (?)
  IMPLICIT NONE
  ! Calling variables
  DOUBLE COMPLEX, INTENT(INOUT) :: workingvol(nx/2+1,ny,nz)
  DOUBLE PRECISION, INTENT(IN)  :: normingvol(nx/2+1,ny,nz)
  INTEGER, INTENT(IN)           :: nx,ny,nz
  ! Other variables
  DOUBLE PRECISION :: norm
  INTEGER :: normterms
  nxby2=nx/2
  ! Obtain the Weiner-filter-like term to the normalizing volume
  normterms=0
  norm=0
  Eachkx: DO kxarray=1,nxby2+1
    Eachky: DO kyarray=1,ny
      Eachkz: DO kzarray=1,nz
        IF (normingvol(kxarray,kyarray,kzarray).ne.0) THEN
          norm=norm+normingvol(kxarray,kyarray,kzarray)
          normterms=normterms+1
        END IF
      END DO Eachkz
    END DO Eachky
  END DO Eachkx
  norm=(0.1*norm)/normterms
  ! Add the Weiner-filter-like term to the normalizing volume
  ! Normalize the working volume
  Wfkx: DO kxarray=1,nx/2+1
    Wfky: DO kyarray=1,ny
      Wfkz: DO kzarray=1,nz
        workingvol(kxarray,kyarray,kzarray)=workingvol(kxarray,kyarray,kzarray)/(normingvol(kxarray,kyarray,kzarray)+norm)
      END DO Wfkz
    END DO Wfky
  END DO Wfkx
  RETURN
  END SUBROUTINE Renormalizemap
  !-------------------------------------------------------
  SUBROUTINE Fsctest(volA,volB,fsc,nx,ny,nz)
  !-------------------------------------------------------
  ! Perform FSC resolution test of two half volumes (volA and volB)
  ! JLR 04/10
  ! JLR 08/10 Fixed bug - only fill fsc(bin) out to bin<=nxby2+1
  ! volA and volB not in wraparound
  IMPLICIT NONE
  ! Calling variables
  INTEGER, INTENT(IN)        :: nx,ny,nz
  DOUBLE COMPLEX, INTENT(IN) :: volA(nx/2+1,ny,nz),volB(nx/2+1,ny,nz)
  REAL, INTENT(OUT)          :: fsc(nx/2+1)
  ! Other
  INTEGER          :: bin
  COMPLEX          :: fsc_num(nx/2+1),term1,term2
  REAL             :: fsc_den1(nx/2+1),fsc_den2(nx/2+1),fsc_den(nx/2+1),radius
  ! The FSC is given by SUM(F1F1*)/SQRT(SUM(|F1|**2)SUM(|F2|**2)) 
  ! The numerator is real because for each F1 that we evaluate there is the Friedel mate F1* that we don't
  ! Therefore: (a+bi)(c+di) + (a-bi)(c-di) = 2ac-2bd for which we just evaluate ac-bd and throw away adi and bci
  ! van Heel, M.  (1987) Ultramicroscopy 21, 95-100.
  nxby2=nx/2
  nyby2=ny/2
  nzby2=nz/2
  ! The following settings to 0 are important!
  fsc_num=0
  fsc_den1=0
  fsc_den2=0
  Eachkz: DO kzarray=1,nz
    Eachky: DO kyarray=1,ny
      Eachkx: DO kxarray=1,nxby2+1
         term1=CMPLX(volA(kxarray,kyarray,kzarray))
         term2=CMPLX(volB(kxarray,kyarray,kzarray))
         kx=REAL(kxarray-1)
         ky=REAL(kyarray-nyby2)
         kz=REAL(kzarray-nzby2)
         radius=SQRT(kx**2+ky**2+kz**2)
         bin=INT(radius)+1
         IF (bin<=nxby2+1) THEN
           IF (ABS(term1)*ABS(term2).ne.0) THEN
             fsc_num(bin) =fsc_num(bin)+term1*CONJG(term2)
             fsc_den1(bin)=fsc_den1(bin)+(ABS(term1))**2
             fsc_den2(bin)=fsc_den2(bin)+(ABS(term2))**2
           END IF
         END IF
      END DO Eachkx
    END DO Eachky
  END DO Eachkz
  Fillfsc: DO kxarray=1,nxby2+1
    fsc_den(kxarray)=SQRT(fsc_den1(kxarray)*fsc_den2(kxarray))
    IF (fsc_den(kxarray).ne.0) THEN
      fsc(kxarray)=REAL(fsc_num(kxarray))/fsc_den(kxarray)
    ELSE
      fsc(kxarray)=0
    END IF
  END DO Fillfsc
  RETURN
  END SUBROUTINE Fsctest
  !--------------------------------------------------------
  REAL FUNCTION Boxft_lut(arg)
  !--------------------------------------------------------
  ! Computes the Fourier transform of the indicator function of a
  ! solid box. The argument ARG should be PI * DSTAR * WIDTH .
  ! This function is modified from frealign
  ! 
  IMPLICIT   NONE
  REAL       sinc, prod, arg(3)
  !         If the argument is large enough, use the standard formula.
  !         If it is small, use the Taylor series expansion.
  prod=1.0
  DO i=1,3
    arg(i)=ABS(arg(i))
    IF (arg(i).LT.2.0E-02) THEN
      sinc=1.0
    ELSE
      j=NINT(180*arg(i))
!      j=NINT(semicirc*arg(i))
      sinc=sinclut(j)
    END IF
    prod=prod*sinc
  END DO
  boxft_lut=prod
  RETURN
  END FUNCTION Boxft_lut

  !-------------------------------------------------------
  SUBROUTINE Moveorigin2d(comparray,nnxby2,nnxby2plus1,ny)
  !-------------------------------------------------------
  IMPLICIT NONE
  ! Calling variables
  DOUBLE COMPLEX, INTENT(INOUT) :: comparray(nnxby2plus1,ny)
  INTEGER, INTENT(IN)           :: nnxby2,nnxby2plus1,ny
  ! Other
  INTEGER        :: ksum
  ! JLR rewrite of optimized code by CM
  Eachky: DO kyarray=1,ny
    ksum=kyarray-1
    IF (IAND(ksum,1).eq.0) THEN
      DO kxarray=2,nnxby2,2
        comparray(kxarray,kyarray)=comparray(kxarray,kyarray)*-1.0
      END DO
    ELSE
      DO kxarray=1,nnxby2plus1,2
        comparray(kxarray,kyarray)=comparray(kxarray,kyarray)*-1.0
      END DO
    END IF
  END DO Eachky
  RETURN
  END SUBROUTINE Moveorigin2d

END MODULE Fouriervolume

MODULE OtherFoMs
  CONTAINS
  !-------------------------------------------------------
  REAL FUNCTION Pres(model,test,shx,shy,nx,ny,rmax1sq,rmax2sq,nxby2,nyby2)
  !------------------------------------------------------
  ! JLR returns an amplitude weighted phase residual between two complex images
  ! Uses the relationships ATAN(X)+ATAN(Y)=ATAN([X+Y]/[1-XY]) and
  !                        ATAX(-X)=-ATAN(X) to calculate a phase difference
  ! Added to fspace_subsandfuncs.f90 16/3/11
  IMPLICIT NONE
  REAL, PARAMETER :: pi = (3.1415926535897)
  INTEGER         :: nxby2,nyby2,kxarray,kyarray,radiussq,lin,kxwrap,kywrap
  INTEGER         :: nx,ny
  REAL            :: rmax1sq,rmax2sq,shx,shy,phase,amp,phasediff,denominator
  DOUBLE COMPLEX  :: model(nxby2+1,ny), test(nxby2+1,ny),modelvox,testvox,pshft,prod
  !
  ! Optimization idea: the sum of amplitudes for a particle does not change
  ! but this quantity is calculated for every particle against every projection
  !
  pres=0
  denominator=0
  shx=2*pi*shx/REAL(nx)
  shy=2*pi*shy/REAL(ny)
  ! Both images are in wraparound order
  DO kyarray=-nyby2+1,nyby2,1 
    DO kxarray=0,nxby2,1
      radiussq=kxarray**2+kyarray**2
      IF(radiussq>=rmax1sq.and.radiussq<=rmax2sq) THEN
        kxwrap=kxarray+1
        lin = kyarray+nyby2                  
        IF (lin<nyby2)  kywrap=lin+nyby2+1
        IF (lin>=nyby2) kywrap=lin-nyby2+1
        phase=shx*REAL(kxarray)+shy*REAL(kyarray)
        pshft=DCMPLX(COS(phase),SIN(phase))
        testvox=test(kxwrap,kywrap)
        amp=ABS(testvox)
        modelvox=(model(kxwrap,kywrap)*pshft)
        prod=testvox*CONJG(modelvox)
        IF(prod.EQ.(0.,0.)) THEN
          phasediff=0.0
        ELSE
          phasediff=ABS(ATAN2(AIMAG(prod),REAL(prod))) ! Uses trig relationships  
        ENDIF
        pres=pres+amp*phasediff
        denominator=denominator+amp
      END IF
    END DO
  END DO
  pres=pres/denominator
  RETURN
  END FUNCTION Pres
END MODULE OtherFoMs
