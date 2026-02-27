 PROGRAM autof2
!------------------------------------------------------------!
! Fortran program used by billy to analyze a list of         !
! energies, reporting which of them is the lowest,           !
! and taking sensible action in the case of                  !
! 2 or more of the energies being the same.                  !
!                                                            !
! Original version : MDT 1992 (second program I ever wrote!) !
! Better version : MDT 2.2002 :-)                            !
!------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER count,num_energies,l,m,n,num,order_of_fit
 DOUBLE PRECISION p(7),e(7),enmin

 open(unit=9,file='dFit.input',status='unknown')
 open(unit=10,file='marker2',status='unknown')
 open(unit=11,file='marker3',status='unknown')


 read(9,*)order_of_fit
 do n=1,7
  read(9,*,err=1,end=1)p(n),e(n)
 enddo
1 enmin=100.d0
 num_energies=n-1
 enmin=minval(e)

 count=0
 do n=1,num_energies
  if(e(n)==enmin)then
   count=count+1
   num=n
  endif
 enddo

 if(count==num_energies)then ! all energies same
  write(10,*)999
  write(11,'(f16.9)')p(1)
  stop
 endif

! If two degenerate energies then pick the value closest to 1 or num_energies.
! (num_energies should be 7 but may not be if something fails to converge).
! If there are two degenerate values equally close to 1 and num_energies, pick
! the first one you find.

 if(count/=0)then

  l=100
  do n=1,num_energies
   if(e(n)==enmin)then
    m=min(abs(n-1),abs(n-num_energies))
    if(m<l)then
     num=n
     l=m
    endif
   endif
  enddo

 endif

 write(10,*)num
 write(11,'(f16.9)')p(num)

 END PROGRAM autof2
