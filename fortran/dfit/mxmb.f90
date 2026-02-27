 SUBROUTINE mxmb(a,nca,nra,b,ncb,nrb,r,ncr,nrr,ncol,nlink,nrow)
!-----------------------------------------------------------------------------!
! General Matrix multiplier (r=r+a*b) taking into account sparsity of b.      !
! r(ncol,nrow)=r(ncol,nrow)+a(ncol,nlink)*b(nlink,nrow)                       !
! Note you need to pre-initialize r.                                          !
!                                                                             !
! V.R. Saunders                                                               !
! Adapted and f90ized MDT.                                                    !
!-----------------------------------------------------------------------------!
 USE dsp
 IMPLICIT NONE
 INTEGER,INTENT(in) :: nca,nra,ncb,nrb,ncr,nrr,ncol,nlink,nrow
 REAL(dp),INTENT(in) :: a(*),b(*)
 REAL(dp),INTENT(inout) :: r(*)
 INTEGER j,ia,ib,ir,iaa,iaa1,iaa2,iaa3,iaa4,iaa5,ibb,irr,iaamax,ncol2,&
  &nca2,ncr2,jaa1,jaa2,jaa3,jaa4,jaa5,loop,jr,krr,kaa
 REAL(dp) fac1,fac2,fac3,fac4,fac5,s1,s2

 ncol2=ncol/2
 iaamax=nlink*nra+1
 nca2=nca+nca ; ncr2=ncr+ncr
 ir=1 ; ib=1

 if((ncol2+ncol2)==ncol)then ! even matrix (number of bands)

  do j=1,nrow

    ibb=ib ; ia=1
1   if(ia==iaamax)goto 15
    fac1=b(ibb)
    iaa1=ia ; ibb=ibb+ncb ; ia=ia+nra
    if(fac1==0.d0)goto 1
    jr=ir+ncr ; jaa1=iaa1+nca
3   if(ia==iaamax)goto 14
    fac2=b(ibb)
    iaa2=ia ; ibb=ibb+ncb ; ia=ia+nra
    if(fac2==0.d0)goto 3
    jaa2=iaa2+nca
5   if(ia==iaamax)goto 13
    fac3=b(ibb)
    iaa3=ia ; ibb=ibb+ncb ; ia=ia+nra
    if(fac3==0.d0)goto 5
    jaa3=iaa3+nca
7   if(ia==iaamax)goto 12
    fac4=b(ibb)
    iaa4=ia ; ibb=ibb+ncb ; ia=ia+nra
    if(fac4==0.d0)goto 7
    jaa4=iaa4+nca
9   if(ia==iaamax)goto 11
    fac5=b(ibb)
    iaa5=ia ; ibb=ibb+ncb ; ia=ia+nra
    if(fac5==0.d0)goto 9
    jaa5=iaa5+nca
    irr=0 ; iaa=0
    do loop=1,ncol2
     s1=r(ir+irr)+fac1*a(iaa1+iaa)+fac2*a(iaa2+iaa)+fac3*a(iaa3+iaa)+ &
      &fac4*a(iaa4+iaa)+fac5*a(iaa5+iaa)
     s2=r(jr+irr)+fac1*a(jaa1+iaa)+fac2*a(jaa2+iaa)+fac3*a(jaa3+iaa)+ &
      &fac4*a(jaa4+iaa)+fac5*a(jaa5+iaa)
     r(ir+irr)=s1
     r(jr+irr)=s2
     irr=irr+ncr2 ; iaa=iaa+nca2
    enddo
    goto 1
11  irr=0 ; iaa=0
    do loop=1,ncol2
     s1=r(ir+irr)+fac1*a(iaa1+iaa)+fac2*a(iaa2+iaa)+fac3*a(iaa3+iaa)+ &
      &fac4*a(iaa4+iaa)
     s2=r(jr+irr)+fac1*a(jaa1+iaa)+fac2*a(jaa2+iaa)+fac3*a(jaa3+iaa)+ &
      &fac4*a(jaa4+iaa)
     r(ir+irr)=s1
     r(jr+irr)=s2
     irr=irr+ncr2 ; iaa=iaa+nca2
    enddo
    goto 15
12  irr=0 ; iaa=0
    do loop=1,ncol2
     s1=r(ir+irr)+fac1*a(iaa1+iaa)+fac2*a(iaa2+iaa)+fac3*a(iaa3+iaa)
     s2=r(jr+irr)+fac1*a(jaa1+iaa)+fac2*a(jaa2+iaa)+fac3*a(jaa3+iaa)
     r(ir+irr)=s1
     r(jr+irr)=s2
     irr=irr+ncr2 ; iaa=iaa+nca2
    enddo
    goto 15
13  irr=0 ; iaa=0
    do loop=1,ncol2
     s1=r(ir+irr)+fac1*a(iaa1+iaa)+fac2*a(iaa2+iaa)
     s2=r(jr+irr)+fac1*a(jaa1+iaa)+fac2*a(jaa2+iaa)
     r(ir+irr)=s1
     r(jr+irr)=s2
     irr=irr+ncr2 ; iaa=iaa+nca2
    enddo
    goto 15
14  irr=0 ; iaa=0
    do loop=1,ncol2
     s1=r(ir+irr)+fac1*a(iaa1+iaa)
     s2=r(jr+irr)+fac1*a(jaa1+iaa)
     r(ir+irr)=s1
     r(jr+irr)=s2
     irr=irr+ncr2 ; iaa=iaa+nca2
    enddo
15  ir=ir+nrr ; ib=ib+nrb

  enddo ! j=1,nrow

 else ! odd matrix (number of bands)

  krr=ncol2*ncr2 ; kaa=ncol2*nca2

  do j=1,nrow

    ibb=ib ; ia=1
16  if(ia==iaamax)goto 29
    fac1=b(ibb)
    iaa1=ia ; ibb=ibb+ncb ; ia=ia+nra
    if(fac1==0.d0)goto 16
17  jr=ir+ncr ; jaa1=iaa1+nca
    if(ia==iaamax)goto 28
    fac2=b(ibb)
    iaa2=ia ; ibb=ibb+ncb ; ia=ia+nra
    if(fac2==0.d0)goto 17
    jaa2=iaa2+nca
19  if(ia==iaamax)goto 27
    fac3=b(ibb)
    iaa3=ia ; ibb=ibb+ncb ; ia=ia+nra
    if(fac3==0.d0)goto 19
    jaa3=iaa3+nca
21  if(ia==iaamax)goto 26
    fac4=b(ibb)
    iaa4=ia ; ibb=ibb+ncb ; ia=ia+nra
    if(fac4==0.d0)goto 21
    jaa4=iaa4+nca
23  if(ia==iaamax)goto 25
    fac5=b(ibb)
    iaa5=ia ; ibb=ibb+ncb ; ia=ia+nra
    if(fac5==0.d0)goto 23
    jaa5=iaa5+nca ; irr=0 ; iaa=0
    do loop=1,ncol2
     s1=r(ir+irr)+fac1*a(iaa1+iaa)+fac2*a(iaa2+iaa)+fac3*a(iaa3+iaa)+ &
      &fac4*a(iaa4+iaa)+fac5*a(iaa5+iaa)
     s2=r(jr+irr)+fac1*a(jaa1+iaa)+fac2*a(jaa2+iaa)+fac3*a(jaa3+iaa)+ &
      &fac4*a(jaa4+iaa)+fac5*a(jaa5+iaa)
     r(ir+irr)=s1
     r(jr+irr)=s2
     irr=irr+ncr2 ; iaa=iaa+nca2
    enddo
    r(ir+krr)=r(ir+krr)+fac1*a(iaa1+kaa)+fac2*a(iaa2+kaa)+fac3*a(iaa3+kaa) &
     &+fac4*a(iaa4+kaa)+fac5*a(iaa5+kaa)
    goto 16
25  irr=0
    iaa=0
    do loop=1,ncol2
     s1=r(ir+irr)+fac1*a(iaa1+iaa)+fac2*a(iaa2+iaa)+fac3*a(iaa3+iaa)+ &
      &fac4*a(iaa4+iaa)
     s2=r(jr+irr)+fac1*a(jaa1+iaa)+fac2*a(jaa2+iaa)+fac3*a(jaa3+iaa)+ &
      &fac4*a(jaa4+iaa)
     r(ir+irr)=s1
     r(jr+irr)=s2
     irr=irr+ncr2 ; iaa=iaa+nca2
    enddo
    r(ir+krr)=r(ir+krr)+fac1*a(iaa1+kaa)+fac2*a(iaa2+kaa)+fac3*a(iaa3+kaa) &
     &+fac4*a(iaa4+kaa)
    goto 29
26  irr=0 ; iaa=0
    do loop=1,ncol2
     s1=r(ir+irr)+fac1*a(iaa1+iaa)+fac2*a(iaa2+iaa)+fac3*a(iaa3+iaa)
     s2=r(jr+irr)+fac1*a(jaa1+iaa)+fac2*a(jaa2+iaa)+fac3*a(jaa3+iaa)
     r(ir+irr)=s1
     r(jr+irr)=s2
     irr=irr+ncr2 ; iaa=iaa+nca2
    enddo
    r(ir+krr)=r(ir+krr)+fac1*a(iaa1+kaa)+fac2*a(iaa2+kaa)+fac3*a(iaa3+kaa)
    goto 29
27  irr=0 ; iaa=0
    do loop=1,ncol2
     s1=r(ir+irr)+fac1*a(iaa1+iaa)+fac2*a(iaa2+iaa)
     s2=r(jr+irr)+fac1*a(jaa1+iaa)+fac2*a(jaa2+iaa)
     r(ir+irr)=s1
     r(jr+irr)=s2
     irr=irr+ncr2 ; iaa=iaa+nca2
    enddo
    r(ir+krr)=r(ir+krr)+fac1*a(iaa1+kaa)+fac2*a(iaa2+kaa)
    goto 29
28  irr=0 ; iaa=0
    do loop=1,ncol2
     s1=r(ir+irr)+fac1*a(iaa1+iaa)
     s2=r(jr+irr)+fac1*a(jaa1+iaa)
     r(ir+irr)=s1
     r(jr+irr)=s2
     irr=irr+ncr2 ; iaa=iaa+nca2
    enddo
    r(ir+krr)=r(ir+krr)+fac1*a(iaa1+kaa)
29  ir=ir+nrr ; ib=ib+nrb

  enddo ! j=1,nrow (2nd time)

 endif ! odd or even matrix

 END SUBROUTINE mxmb
