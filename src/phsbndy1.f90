SUBROUTINE phsbndy1(px,hw,hs)
  !     Purpose:  To calculate the enthaply of sat'd water and sat'd steam
  !             for a given pressure.
  !     Note: phsbndy1 is identical to phsbndy2 except that fewer parameters
  !             are calculated.
  !     Method: This subroutine uses a second derivative cubic spine.
  !     Reference:  Press et al., 1986, Numerical Recipes, the Art of
  !                   Scientific Computing, Chapter 3.3 on Cubic Spline
  !                   Interpolation, p. 86-87.
  !     Range: This spline covers the range in pressure from 0.5 bars
  !               to the critical point, 220.55 bars with 134 selected
  !               pressure-enthalpy pairs.  The pressure enthalpy pairs
  !               were calculated from the steam tables (Haar et al., 1984)
  USE machine_constants, ONLY: kdp
  USE tables
  IMPLICIT NONE
  REAL(kind=kdp), INTENT(IN) :: px
  REAL(kind=kdp), INTENT(OUT) :: hw, hs
  !       PX - pressure in dyne/cm^2 
  !       HW - enthalpy of saturated water in erg/g
  !       HS - enthalpy of saturated steam in erg/g
  INTEGER :: khi, klo
  REAL(kind=kdp) :: a, a3, b, b3, h, h2
  REAL(kind=kdp) :: pxx
  !
  !       Ps - saturation pressures at selected points
  !       Swh - saturated water enthalpy at selected points
  !       D2swh - saturated water enthalpy second derivative w.r.t. pressure
  !       Ssh - saturated steam enthalpy at selected points
  !       D2ssh - saturated steam enthalpy second derivative w.r.t. pressure
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3 $//$Date: 2001/08/10 19:01:51 $'
  !     ------------------------------------------------------------------
  !...
  !.....Special patch for low pressures
  pxx=MAX(px,5.e5_kdp)
  CALL phsbnklo (1, pxx, klo)
  !---check to be sure that klo pointer is positioned correctly
  DO WHILE (pxx > ps(klo+1) .AND. klo < mdimws-1)
     klo = klo + 1
  END DO
  DO WHILE (pxx < ps(klo) .AND. klo > 1)
     klo = klo - 1
  END DO
  !---Interpolate with the cubic spline
  khi = klo + 1
  h = ps(khi) - ps(klo)
  a = (ps(khi)-pxx)/h
  b = (pxx-ps(klo))/h
  a3 = a**3 - a
  b3 = b**3 - b
  h2 = h**2/6.0_kdp
  hw = a*swh(klo) + b*swh(khi) + (a3*d2swh(klo)+b3*d2swh(khi))*h2
  hs = a*ssh(klo) + b*ssh(khi) + (a3*d2ssh(klo)+b3*d2ssh(khi))*h2
END SUBROUTINE phsbndy1
