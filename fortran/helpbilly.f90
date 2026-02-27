 PROGRAM autof
!------------------------------------------------------------!
! Fortran routine used by billy to compare real numbers      !
!                                                            !
! MDT 1992 - first fortran program I ever wrote!             !
!------------------------------------------------------------!
 IMPLICIT NONE
 DOUBLE PRECISION emin,emin0
 open(unit=9,file='eTemp',status='unknown')
 open(unit=10,file='marker',status='unknown')
 read(9,*)emin
 read(9,*)emin0
 if(emin<emin0)then
  write(10,*)"1"
 else
  write(10,*)"0"
 endif
END PROGRAM autof
