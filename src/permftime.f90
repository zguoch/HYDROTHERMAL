SUBROUTINE permftime
  !     Purpose:  To compute the permeability
  !               of selected rock types as a function of time.
  !
  USE machine_constants, ONLY: kdp
  USE mesh
  USE parameters
  IMPLICIT NONE
  INTEGER :: ic, icdum, irx, irxok
  INTEGER :: ndim
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.1 $//$Date: 2002/03/18 23:38:08 $'
  !     ------------------------------------------------------------------
  ! ...    determine the rock type number for the first active node
  irxok = icrxtype(npic(1))
  ! ...        ------------ X Permeability ---------------
  ndim = nxyz
  IF (irxftimopt(2,irxok,1) > 0) CALL permfunc(xk,ndim,2)
  ! ...                 ----- Y Permeability -----
  IF (irxftimopt(3,irxok,1) > 0) THEN
     CALL permfunc(yk,ndim,3)
  ELSE IF (irxftimopt(2,irxok,1) > 0 .AND. ykrxkxfac(1) /= 0._kdp) THEN
     DO icdum = 1,npiccons
        ic = npic(icdum)
        irx = icrxtype(ic)
        yk(ic) = ykrxkxfac(irx)*xk(ic)
     END DO
  END IF
  ! ...                 ----- Z Permeability -----
  IF (irxftimopt(4,irxok,1) > 0) THEN
     CALL permfunc(zk,ndim,4)
  ELSE IF (irxftimopt(2,irxok,1) > 0 .AND. zkrxkxfac(1) /= 0._kdp) THEN
     DO icdum = 1,npiccons
        ic = npic(icdum)
        irx = icrxtype(ic)
        zk(ic) = zkrxkxfac(irx)*xk(ic)
     END DO
  END IF
END SUBROUTINE permftime
