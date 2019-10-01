SUBROUTINE resid
  ! ... Purpose:  To compute the residual mass and energy, F equations, 
  ! ...       for each active node.
  ! ...       Also compute the Linf-norm of Fp and Fh
  USE machine_constants, ONLY: kdp
  USE bc
  USE fdeq
  USE mesh
  USE source
  USE variables
  IMPLICIT NONE
  INTEGER :: i, ic, ik, j, k
  REAL(KIND=kdp) :: vol
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3.1.5 $//$Date: 2007/08/17 18:29:39 $'
  !     ------------------------------------------------------------------
  resmaxm = 0._kdp
  resmaxe = 0._kdp
  residm = 0._kdp
  reside = 0._kdp
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     ! ... Skip if node is a constant P and H cell or is not in active domain
     ! ...    or is a seeping cell
     ! ...    Energy residual for seeping cell is missed****
     IF(ibc(ik,j) == -1 .OR. ibc(ik,j) == 11 .OR.  &
          ((ibc(ik,j) == 30 .OR. ibc(ik,j) == 50) .AND. seeping(ik,j))) CYCLE
     ! ... Calculate spatial flow terms of governing finite difference equations
     CALL goveqn(residm(ic),reside(ic),i,j,k)
     ! ... add in the capacitance and source terms
     residm(ic) = residm(ic) - (xm(ic)-xmoldt(ic))/delt + q(ic)
     reside(ic) = reside(ic) - (en(ic)-enoldt(ic))/delt + qh(ic)
     ! ... Cancel residual for specified pressure nodes
     IF (ibc(ik,j) == 10) residm(ic) = 0._kdp
     ! ... Convert to specific residual (per unit volume)
     vol = dx(i)*dy(j)*dz(k)
     IF(irad) vol = 3.1415927_kdp*(rsq(i+1)-rsq(i))*dz(k)
     residm(ic) = residm(ic)/vol
     reside(ic) = reside(ic)/vol
     ! ... Store value and location of maximum residuals
     ! ... Calculate the Linf-norm of Fp and Fh; with sign
     IF (ABS(residm(ic)) > ABS(resmaxm)) THEN
        resmaxm = residm(ic)
        imres = i
        jmres = j
        kmres = k
     END IF
     IF (ABS(reside(ic)) > ABS(resmaxe)) THEN
        resmaxe = reside(ic)
        ieres = i
        jeres = j
        keres = k
     END IF
  END DO
END SUBROUTINE resid

