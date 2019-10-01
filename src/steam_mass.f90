SUBROUTINE steam_mass(wvap, wliq, wscrit)
  ! ... Purpose: To determine the mass of steam in the region
  ! ...          and the mass of sup
  USE machine_constants, ONLY: kdp
  USE math_constants
  USE bc, ONLY: ibc
  USE mesh
  USE parameters
  USE variables
  IMPLICIT NONE
  REAL(KIND=kdp), INTENT(OUT) :: wvap, wliq, wscrit
  !
  INTEGER :: i, ic, ik, j, k
  REAL(KIND=kdp) :: vol
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4 $//$Date: 2006/10/03 23:35:36 $'
  !     ------------------------------------------------------------------
  wvap = 0._kdp
  wliq = 0._kdp
  wscrit = 0._kdp
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF (ibc(ik,j) /= -1) THEN        ! ... Select active cells
        vol = dx(i)*dz(k)*dy(j)
        IF (irad) vol = pi*(rsq(i+1)-rsq(i))*dz(k)
        IF (p(ic) < 220.55e6_kdp .AND. ind(ic) /= 5)  THEN
           wvap = wvap + dens(ic)*(1._kdp-satnw(ic))*vol*phi(ic)
           wliq = wliq + denw(ic)*satnw(ic)*vol*phi(ic)
        ELSEIF(p(ic) >= 220.55e6_kdp) THEN
           wscrit = wscrit + denw(ic)*vol*phi(ic)     ! ... denw=dens in supercritical region
        ELSEIF(ind(ic) == 5) THEN
           wliq = wliq + denw(ic)*satnw(ic)*vol*phi(ic)
        END IF
     END IF
  END DO
END SUBROUTINE steam_mass
