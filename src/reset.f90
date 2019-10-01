SUBROUTINE reset
  ! ... Purpose:  To reset all variables to the values they had at
  ! ...           the end of the previous time step.  This
  ! ...           subroutine is only called if the previous sequence of
  ! ...           N-R iterations did not converge or the maximum allowed
  ! ...           change in saturation was exceeded.  A smaller time step
  ! ...           can then be used to attempt a convergent sequence of N-R iterations.
  USE machine_constants, ONLY: kdp
  USE bc
  USE i_c
  USE fdeq, ONLY: ioptupst
  USE mesh
  USE parameters
  USE variables
  IMPLICIT NONE
  INTEGER :: i, ic, ik, j, k, lll
  INTERFACE
     SUBROUTINE tempdegc(p,h,tc)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
       REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: tc
     END SUBROUTINE tempdegc
  END INTERFACE
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.9 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  DO  lll=1,npiccons
     ic = npic(lll)
     p(ic) = poldt(ic)
     h(ic) = holdt(ic)
     ind(ic) = indoldt(ic)
     indoldnr(ic) = indoldt(ic)
     xm(ic) = xmoldt(ic)
     en(ic) = enoldt(ic)
     ! ...   Xmoldt(ic) = Xmoldt2(ic)
     ! ...   Enoldt(ic) = Enoldt2(ic)
     phi(ic) = phiinit(ic)*(1._kdp+(p(ic)-pinit(ic))*beta(ic))
  END DO
  ! ... Now do the reset of pressure values preserving the seeping cell
  ! ...      atmospheric pressure values
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     ! ... Select seepage nodes that are seeping
     IF ((ibc(ik,j) == 30 .OR. ibc(ik,j) == 50) .AND. seeping(ik,j)) THEN
        p(ic) = patm
     ELSE
        p(ic) = poldt(ic)
     END IF
  END DO
  ! ... Recalculate temperature using reset values of P and H
  ! ...        for permeability calculations and rock enthalpy calculations
  CALL tempdegc(p,h,tc)
  IF (lrxftreq) THEN       ! ... Calculate permeability and transmissivity as f(Temperature)
     CALL rockparm
     CALL tcalc     ! ... Calculate transmissivity and thermal conductivity for cell faces
  END IF
  IF (lrxftimreq) THEN       ! ... Calculate permeability and transmissivity as f(Time) 
     CALL permftime
     CALL tcalc     ! ... Calculate transmissivity and thermal conductivity for cell faces
  END IF
  CALL enthrock      ! ... Calculate enthalpy of the porous matrix
  CALL prpty(1)      ! ... Calculate fluid property coefficients and derivatives
  CALL wellallo      ! ... Allocate the heat and mass to the source/sink terms
  IF (ioptupst) CALL upstre  ! ... Determine the upstream node for each cell face
END SUBROUTINE reset
