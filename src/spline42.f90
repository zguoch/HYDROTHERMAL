SUBROUTINE spline42(h1, d1, pd1ph, pd1pp, t1, pt1ph, pt1pp, v1,  &
     pv1ph, pv1pp, h2, d2, pd2ph, pd2pp, t2, pt2ph,  &
     pt2pp, v2, pv2ph, pv2pp, hx, da, pdaph, pdapp,  &
     ta, ptaph, ptapp, va, pvaph, pvapp)
  !     Purpose:  to interpolate the density, temperature, and viscosity of
  !          water/steam, and the derivatives of those values w.r.t. enthalpy
  !          and pressure at a specific enthalpy-pressure point between
  !          two points with constant pressure but different enthalpies.

  !     Method: the density, temperature and viscosity values are interpolated
  !          using a modified cubic spline in which the values and the
  !          the derivatives w.r.t. enthalpy at the two end points are used.
  !          The modified cubic spline for two points is a combination of
  !          the subroutines SPLINE and SPLINT from Press et al., 1987,
  !          Numerical Recipies, p. 88-89.
  !          A simple linear interpolation is used to determine the derivatives
  !          w.r.t. enthalpy and pressure.
  !     NOTE:  H1 must be > H2
  USE machine_constants, ONLY: kdp
  IMPLICIT NONE
  REAL(kind=kdp), INTENT(IN) :: h1, d1, pd1ph, pd1pp, t1, pt1ph, pt1pp, v1, pv1ph,  &
       pv1pp, h2, d2, pd2ph, pd2pp, t2, pt2ph, pt2pp, v2, pv2ph, pv2pp, hx
  REAL(kind=kdp), INTENT(OUT) :: da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp
  !     Passed variables:
  !       H1 - (input) enthalpy at H,P point 1 - in ergs/g (1 J = 1.0D7 ergs)
  !       H2 - (input) enthalpy at H,P point 2 - in ergs/g (1 J = 1.0D7 ergs)
  !       HX - (input) enthalpy at H,P point X - in ergs/g (1 J = 1.0D7 ergs)
  !       D# - fluid density at H,P point # (g/cm**3)
  !       PD#PH - derivative of density w.r.t. enthalpy at H,P point #
  !       PD#PP - derivative of density w.r.t. pressure at H,P point #
  !       T# - temperature at H,P point # (deg C)
  !       PT#PH - derivative of temperature w.r.t. enthalpy at H,P point #
  !       PT#PP - derivative of temperature w.r.t. pressure at H,P point #
  !       V# - viscosity at H,P point # ()
  !       PV#PH - derivative of viscosity w.r.t. enthalpy at H,P point #
  !       PV#PP - derivative of viscosity w.r.t. pressure at H,P point #
  !     where # = 1,2,a
  !         1 and 2 -  refer to the respective input values at H1 and H2
  !         a  - refers to the returned values interpolated for point HX
  !
  REAL(kind=kdp) :: a, b, h, u1d, u1t, u1v, u2d, u2t, u2v,  &
       z1d, z1t, z1v, z2d, z2t, z2v
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.2 $//$Date: 2000/12/12 21:27:01 $'
  !     ------------------------------------------------------------------
  h = h1 - h2
  a = (h1-hx)/h
  b = (hx-h2)/h
  !---density
  u2d = (3._KDP/h)*((d1-d2)/h-pd2ph)
  u1d = (3._KDP/h)*(pd1ph-(d1-d2)/h)
  z1d = (u1d-0.5_KDP*u2d)/0.75D0
  z2d = -0.5D0*z1d + u2d
  da = a*d2 + b*d1 + ((a**3-a)*z2d+(b**3-b)*z1d)*(h**2) /6._KDP
  pdaph = a*pd2ph + b*pd1ph
  pdapp = a*pd2pp + b*pd1pp
  !---temperature
  u2t = (3._KDP/h)*((t1-t2)/h-pt2ph)
  u1t = (3._KDP/h)*(pt1ph-(t1-t2)/h)
  z1t = (u1t-0.5_KDP*u2t)/0.75D0
  z2t = -0.5D0*z1t + u2t
  ta = a*t2 + b*t1 + ((a**3-a)*z2t+(b**3-b)*z1t)*(h**2) /6._KDP
  ptaph = a*pt2ph + b*pt1ph
  ptapp = a*pt2pp + b*pt1pp
  !---viscosity
  u2v = (3._KDP/h)*((v1-v2)/h-pv2ph)
  u1v = (3._KDP/h)*(pv1ph-(v1-v2)/h)
  z1v = (u1v-0.5_KDP*u2v)/0.75_kdp
  z2v = -0.5_kdp*z1v + u2v
  va = a*v2 + b*v1 + ((a**3-a)*z2v+(b**3-b)*z1v)*(h**2) /6._KDP
  pvaph = a*pv2ph + b*pv1ph
  pvapp = a*pv2pp + b*pv1pp
END SUBROUTINE spline42
