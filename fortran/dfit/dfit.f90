  PROGRAM dfit
!-----------------------------------------------------!
!   DFIT - Simple fitting of polynomial to data       !
!   V.R.Saunders                                      !
!   Translated to F90 and formatted, MDT.             !
!-----------------------------------------------------!
  USE dsp
  IMPLICIT NONE
  INTEGER,PARAMETER :: nmax=20
  INTEGER i,j,k,n,iord(20),is(nmax),nord,nord1,ierr
  REAL(dp) x(20),y(20),a(20,21),s(nmax),x0,y0,x1,y1,eps,det,val,vrs,al,fac
  INTERFACE
   SUBROUTINE lesolv(a,b,work,ipivot,id,n,nrhs)
    USE dsp
    INTEGER n,id,nrhs,ipivot(n)
    REAL(dp) a(id,n),b(id,nrhs),work(n)
   END SUBROUTINE lesolv
  END INTERFACE

  read(5,*)nord
  do i=1,nord
   iord(i)=i-1
  enddo
  write(6,'(////'' LENGTH OF POLYNOMIAL='',i4//'' POWERS IN POLYNOMIAL='', &
   &20i3)')nord,(iord(i),i=1,nord)

  n=0
  do
   n=n+1
   read(5,*,iostat=ierr)x(n),y(n)
   if(ierr<0)exit
   write(6,'('' X,Y='',f17.9,e21.12)')x(n),y(n)
  enddo

  n=n-1
  nord1=nord+1
  x0=x(1)
  y0=y(1)

  a(1:nord,1:nord1)=0.d0

  do i=1,n
   x1=x(i)-x0
   y1=y(i)-y0
   do j=1,nord
    s(j)=1.d0
    if(iord(j)/=0)s(j)=x1**iord(j)
    a(j,nord1)=s(j)*y1+a(j,nord1)
    do k=1,j
     a(k,j)=s(j)*s(k)+a(k,j)
     a(j,k)=a(k,j)
    enddo ! k
   enddo ! j
  enddo ! i

  eps=1.d-40
  call lesolv(a,a(1,nord1),s,is,nmax,nord,1)
  write(6,'('' COEFFICIENTS='',(/1e21.12))')(a(i,nord1),i=1,nord)

  det=0.d0
  do i=1,n
   x1=x(i)-x0
   val=0.d0
   do j=1,nord
    vrs=1.d0
    if(iord(j)/=0)vrs=x1**iord(j)
    val=val+vrs*a(j,nord1)
   enddo ! j
   val=val+y0
   det=(y(i)-val)**2+det
   write(6,'('' X,Y,YFIT='',f16.9,2e21.12)')x(i),y(i),val
  enddo ! i

  det=sqrt(det/dble(n))
  write(6,'('' RMS DEVIATION='',1pe15.5)')det

  al=0.d0
  do k=1,9
   vrs=0.d0
   det=0.d0
   fac=1.d0
   do i=2,nord
    s(i-1)=a(i,nord1)*(i-1)
    vrs=fac*s(i-1)+vrs
    fac=fac*al
   enddo ! i
   fac=1.d0
   do i=3,nord
    det=s(i-1)*fac*(i-2)+det
    fac=fac*al
   enddo
   if(det/=0.d0)then
    al=al-vrs/det
   else
    write(6,*)'Division by zero.'
    stop
   endif
  enddo ! k

  vrs=al+x0
  write(6,'(/'' MINIMUM AT X='',f16.7)')vrs
  vrs=y0+a(1,nord1)
  det=al
  do i=2,nord
   vrs=det*a(i,nord1)+vrs
   det=det*al
  enddo
  write(6,'(/'' Y AT MINIMUM='',e21.12)')vrs

  END PROGRAM dfit
