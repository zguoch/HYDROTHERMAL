SUBROUTINE relperm2 (px, hx, xkw, pxkwp, pxkwh, xks, pxksp, pxksh)
  !     Purpose:  To calculate the relative permeability of "supercritical water"
  !               in order to smooth the transition from water and steam
  !               via the supercritical route
  !
  !     Method:  The relative permeabilities vary from 1.0 on the water
  !              side of the critical point to 0.0 on the steam side.
  !              Relative permeabilities decrease in an arcuate manner
  !              around a point (hub) located just inside the two-phase
  !              region at P= 220.0 bars and H= 2086 J/g
  USE machine_constants, ONLY: kdp
  IMPLICIT NONE
  REAL(kind=kdp), INTENT(IN) :: px, hx
  REAL(kind=kdp), INTENT(OUT) :: xkw, pxkwp, pxkwh, xks, pxksp, pxksh
  !     Passed Variables:
  !       PX (input) - pressure (dyne/cm^2)
  !       HX (input) - enthalpy (erg/g)
  !       XKW (retn) - relative permeability of water
  !       PXKWP (retn) - deriv. of relative permeability of water w.r.t. P
  !       PXKWH (retn) - deriv. of relative permeability of water w.r.t. H
  !       XKS (retn) - relative permeability of steam
  !       PXKSP (retn) - deriv. of relative permeability of steam w.r.t. P
  !       PXKSH (retn) - deriv. of relative permeability of steam w.r.t. H
  !
  REAL(kind=kdp) :: aaa, daaah, daaap, hhub, phub
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.2 $//$Date: 2000/12/12 21:27:01 $'
  !     ------------------------------------------------------------------
  phub = 220.0e6_kdp
  hhub = 20.86e9_kdp
  !---if pressure is less than the "hub" pressure
  IF (px <= phub) THEN
     IF (hx < hhub) xkw = 1._KDP
     IF (hx > hhub) xkw = 0._KDP
     pxkwp = 0._KDP
     pxkwh = 0._KDP
     !---if pressure is greater than the "hub" pressure
     !---if enthalpy equals "hub" enthalpy
  ELSE IF (hx == hhub) THEN
     xkw = 0.5_kdp
     pxkwp = 0._kdp
     pxkwh = 0._kdp
  ELSE
     aaa = (px-phub)/(hx-hhub)
     daaap = 1._KDP/(hx-hhub)
     daaah = -(px-phub)/((hx-hhub)*(hx-hhub))
     !---if enthalpy is less than the  "hub" enthalpy
     IF (hx < hhub) THEN
        xkw = (0.5_kdp*COS(ATAN(aaa))) + 0.5_kdp
        pxkwp = -0.5_kdp*SIN(ATAN(aaa))*(1._kdp/(1._kdp+aaa*aaa))*daaap
        pxkwh = -0.5_kdp*SIN(ATAN(aaa))*(1._kdp/(1._kdp+aaa*aaa))*daaah
        !---if enthalpy is greater than the  "hub" enthalpy
     ELSE
        xkw = (-0.5_kdp*COS(ATAN(aaa))) + 0.5_kdp
        IF (xkw < 1.0e-20_kdp) xkw = 0._kdp
        pxkwp = 0.5_kdp*SIN(ATAN(aaa))*(1._kdp/(1._kdp+aaa*aaa))*daaap
        pxkwh = 0.5_kdp*SIN(ATAN(aaa))*(1._kdp/(1._kdp+aaa*aaa))*daaah
     END IF
  END IF
  xks = 1._kdp - xkw
  pxksp = -pxkwp
  pxksh = -pxkwh
END SUBROUTINE relperm2
