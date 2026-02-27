  SUBROUTINE lesolv(a,b,work,ipivot,id,n,nrhs)
!------------------------------------------------------------------!
! Linear equation solver by LU decomposition.                      !
! V.R. Saunders 1991                                               !
! Translated into F90 and formatted, MDT.                          !
!------------------------------------------------------------------!
  USE dsp
  IMPLICIT NONE
  INTEGER i,j,k,n,id,jm1,jp1,nrhs,ipivot(n)
  REAL(dp) a(id,n),b(id,nrhs),work(n),bilbo,frodo
  REAL(dp),PARAMETER :: toll=1.d-30
  INTERFACE
   SUBROUTINE mxmb(a,nca,nra,b,ncb,nrb,r,ncr,nrr,ncol,nlink,nrow)
    USE dsp
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nca,nra,ncb,nrb,ncr,nrr,ncol,nlink,nrow
    REAL(dp),INTENT(in) :: a(*),b(*)
    REAL(dp),INTENT(inout) :: r(*)
   END SUBROUTINE mxmb
   SUBROUTINE errvrs(ierr,namz,mzss)
    INTEGER ierr
    CHARACTER(6) namz
    CHARACTER(63) mzss
   END SUBROUTINE errvrs
  END INTERFACE

! Decomposition step
  do j=1,n
   ipivot(j)=j
  enddo
  do j=1,n
   jm1=j-1
   jp1=j+1
   if(jm1/=0)call mxmb(a(j,1),1,id,a(1,j),1,id,a(j,j),1,id,n-jm1,jm1,1)
   bilbo=abs(a(j,j))
   k=0
   do i=jp1,n
    frodo=abs(a(i,j))
    if(bilbo>=frodo)cycle
    bilbo=frodo
    k=i
   enddo
   if(bilbo<toll)call errvrs(0,'LESOLV','THE MATRIX ON THE LHS OF THE &
    &EQUATION IS SINGULAR              ')
   if(k/=0)then
    i=ipivot(j)
    ipivot(j)=ipivot(k)
    ipivot(k)=i
    do i=1,n
     bilbo=a(j,i)
     a(j,i)=a(k,i)
     a(k,i)=bilbo
    enddo
   endif
   bilbo=-1.d0/a(j,j)
   a(j,j)=bilbo
   if(j/=n)then
    if(jm1/=0)call mxmb(a(1,jp1),id,1,a(j,1),id,1,a(j,jp1),id,1,n-j,jm1,1)
    do i=jp1,n
     a(j,i)=a(j,i)*bilbo
    enddo
   endif
  enddo

  do k=1,nrhs
! Permute the right hand side
   do i=1,n
    work(i)=b(ipivot(i),k)
   enddo
! Forward elimination step
   do j=1,n
    bilbo=a(j,j)*work(j)
    work(j)=-bilbo
    do i=j+1,n
     work(i)=a(i,j)*bilbo+work(i)
    enddo
   enddo
! Back substitution step
   do j=n,2,-1
    bilbo=work(j)
    do i=1,j-1
     work(i)=a(i,j)*bilbo+work(i)
    enddo
   enddo
   do i=1,n
    b(i,k)=work(i)
   enddo
  enddo

  END SUBROUTINE lesolv
