  SUBROUTINE errvrs(ierr,namz,mzss)
  INTEGER ierr,iout
  CHARACTER(6) namz
  CHARACTER(63) mzss
  iout=0
  if(ierr==0)then
   write(iout,1)namz,mzss
1  format(' ERROR **** ',a6,' **** ',a63)
   stop
  else
   write(iout,2)namz,mzss
2  format(' WARNING **** ',a6,' **** ',a63)
  endif
  END SUBROUTINE errvrs
