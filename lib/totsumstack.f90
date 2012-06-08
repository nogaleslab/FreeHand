!************************************************
!------------------------------------------------
! v1.00     JLR 11/02
!-------------------------------------------------
! INPUTS:
! CARD 1: (INFILE)  Stack containing particles MRC Format
! CARD 2: (OUTFILE) Average of all images in stack
!-------------------------------------------------


!*************************************************
!*************************************************
!*************************************************

PROGRAM Totsum
IMPLICIT NONE
INTEGER :: nparts,nx,ny
INTEGER :: nxyzi(3),mxyzi(3),nxyzf(3),modei
INTEGER :: nxyz(3),mxyz(3),mode,nxyzst(3)
REAL    :: dmini,dmaxi,dmeani
REAL    :: dmin,dmax,dmean,cell(6)
REAL,ALLOCATABLE :: stack(:,:,:)              ! input stack
REAL,ALLOCATABLE :: outimage(:,:)             ! output average
INTEGER :: sec,lin,vox
CHARACTER(80) :: infile,outfile
CHARACTER(70) :: title

READ (5,*) infile
READ (5,*) outfile

CALL Imopen (1,infile,"OLD")
CALL IRDHDR(1,nxyzi,mxyzi,modei,dmini,dmaxi,dmeani)

nx     = nxyzi(1)
ny     = nxyzi(2)
nparts = nxyzi(3)
nxyz  = (/nx,ny,1/) 

ALLOCATE(stack(nx,ny,nparts))      ;   stack     = 0
ALLOCATE(outimage(nx,ny))          ;   outimage  = 0

! Put input stack into memory
Eachsection: DO sec=0,nparts-1
  CALL IRDSEC(1,stack(:,:,sec+1),*999)
END DO Eachsection

WRITE (*,*) "Input stack is:",infile
WRITE(*,*)  "Stack put into array"

DO sec = 1,nparts
  outimage  = outimage + stack(:,:,sec)/nparts  
END DO

CALL Imopen(2,outfile,"UNKNOWN")
dmin=MINVAL(outimage)
dmax=MAXVAL(outimage)
dmean=SUM(outimage)/(nx*ny)


nxyz=(/nx,ny,1/)
nxyzst=(/0,0,0/)
mxyz=(/nx,ny,1/)
cell=(/nx,ny,1,90,90,90/)



CALL Icrhdr(2,nxyz,mxyz,2,"   ",1)
CALL Ialsiz(2,nxyz,nxyzst)
CALL Itrcel(2,1)
CALL Ialcel(2,cell)
CALL Irtcel(2,cell)
CALL Iwrhdr(2,title,-1,dmin,dmax,dmean)
CALL Iwrsec(2,outimage) !! Write the sum image to the output file
CALL Imclose(2)
CALL Imclose(1)

STOP
999 STOP "End of file read error"
END PROGRAM Totsum
