SUBROUTINE storativ
  ! ... Purpose: Compute the fluid mass and energy stored in each cell
  USE machine_constants, ONLY: kdp
  USE bc, ONLY: ibc
  USE math_constants
  USE mesh
  USE parameters
  USE variables
  IMPLICIT NONE
  INTEGER :: i, ic, ik, j, k
  REAL(KIND=kdp) :: avden, hs, hw, vol
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.8 $//$Date: 2006/10/03 23:37:49 $'
  !     ------------------------------------------------------------------
  xm = 0._kdp
  en = 0._kdp
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF (ibc(ik,j) /= -1) THEN        ! ... Select active cells
        vol = dx(i)*dz(k)*dy(j)
        IF (irad) vol = pi*(rsq(i+1)-rsq(i))*dz(k)
        IF(ind(ic) /= 5) THEN
           avden = denw(ic)*satnw(ic) + dens(ic)*(1._kdp-satnw(ic))
        ELSE IF(ind(ic) == 5) THEN
           avden = denw(ic)*satnw(ic)
        END IF
        ! ... mass in storage
        xm(ic) = vol*phi(ic)*avden
        ! ... energy in storage
        IF (ind(ic) == 2) THEN     ! ... Two phase water-steam
           CALL phsbndy1(p(ic),hw,hs)
           en(ic) = vol*(phi(ic)*(satnw(ic)*denw(ic)*hw  &
                    + (1._kdp-satnw(ic))*dens(ic)*hs)  &
                    + (1._kdp-phi(ic))*hrock(ic)*df(ic))
        ELSE                       ! ... Single phase and air-water two phase
           en(ic) = vol*(avden*phi(ic)*h(ic)  &
                    + (1._kdp-phi(ic))*hrock(ic)*df(ic))
        END IF
     END IF
  END DO
END SUBROUTINE storativ
