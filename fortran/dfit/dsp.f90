MODULE dsp
!---------------------------------------------!
! Constants for variable declarations.        !
!---------------------------------------------!
IMPLICIT NONE

! e.g. Number of bytes for integer type n such that -10^9 < n < 10^9.
! INTEGER,PARAMETER :: i9b=selected_int_kind(9)

! Number of bytes for default single precision type.
INTEGER,PARAMETER :: sp=kind(1.)

! Number of bytes for default double precision type.
INTEGER,PARAMETER :: dp=kind(1.d0)

END MODULE dsp
