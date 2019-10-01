FUNCTION interp1(xarg,nx,xs,fs) RESULT(xinterp)
  ! ... Linear 1-d interpolation
  ! ... XS must be in algebraic ascending order
  ! ... If XARG is outside the range of data, extrapolation is done
  ! ...      for one dimensional invocation
  USE machine_constants, ONLY: kdp
  IMPLICIT NONE
  REAL(KIND=kdp), INTENT(IN) :: xarg
  INTEGER, INTENT(IN) :: nx
  REAL(KIND=kdp), DIMENSION(nx), INTENT(IN) :: xs, fs
  REAL(KIND=kdp) :: xinterp
  !
  INTEGER :: i
  REAL(KIND=kdp) :: a1
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3 $//$Date: 2007/10/16 21:17:22 $'
  !     ------------------------------------------------------------------
  !...
  ! ... Find location of xarg in table
  DO i=1,nx
     IF(xs(i) > xarg) EXIT
  END DO
  i = max(i,2)
  i = min(i,nx)
  a1 = (xarg-xs(i-1))/(xs(i)-xs(i-1))
  xinterp = (1._kdp-a1)*fs(i-1) + a1*fs(i)
END FUNCTION interp1
