SUBROUTINE bndyt0(px, hx, t0d1, pt0dph1, pt0dpp1, t0t1, pt0tph1,  &
     pt0tpp1, t0v1, pt0vph1, pt0vpp1)
  ! ... Purpose:  To calculate the enthalpy, temperature, density, and
  ! ...     viscosity of water at a temperature of 0 Deg.C for a given
  ! ...     pressure (greater than 240 bars).  Also
  ! ...     calculates the derivatives of T,D, and V w.r.t. H and P.
  ! ...     A cubic spline is used to determine the values.
  ! ...     This spline covers the range in pressure from 240 - 10,000 bars
  ! ...     with 100 selected pressure-enthalpy pairs
  ! ...     Reference:  Press et al., 1986, Numerical Recipes, the Art of
  ! ...          Scientific Computing, Chapter 3.3 on Cubic Spline
  ! ...          Interpolation, p. 86-87.
  USE machine_constants, ONLY: kdp
  USE f_units
  USE control
  USE tables
  IMPLICIT NONE
  REAL(KIND=kdp), INTENT(IN) :: px
  REAL(KIND=kdp), INTENT(OUT) :: hx, t0d1, pt0dph1, pt0dpp1, t0t1, pt0tph1, pt0tpp1, &
       t0v1, pt0vph1, pt0vpp1
  ! ...   PX - pressure [dynes/cm^2] 
  ! ...   HX - enthalpy of saturated steam [ergs/g] 
  !
  INTEGER :: khi, klo
  REAL(KIND=kdp) :: a, a3, b, b3, h, h2, logp
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5 $//$Date: 2007/07/31 19:52:56 $'
  !     ------------------------------------------------------------------
  !...
  ! ... if pressure is out of range of cubic spline, abort
  IF (px < 240.e6_kdp .OR. px > 1.0e10_kdp) THEN
     WRITE(fustdout, '(A/10X,A,1P,G13.6)')  &
          ' **** ERROR:Pressure beyond limits of cubic spline in bndyT0' , 'P= ', px
     WRITE(fuclog, '(A/10X,A,1P,G13.6)')  &
          ' **** ERROR:Pressure beyond limits of cubic spline in bndyT0' , 'P= ', px
     ierr(70) = .TRUE.
     RETURN
  END IF
  ! ... estimate pressure position in the lookup table
  logp = LOG10(px)
  klo = INT(logp*50._kdp) - 418  ! ... for 240 <= P <~ 7,000 bars step by 0.02 log bars
  IF (logp >= 9.84_kdp .AND. logp < 9.9_kdp) THEN
     klo = INT(logp*100._kdp) - 910  ! ... for 7,000 < P <= 8,000 bars - step by 0.01 log bars
  ELSEIF (logp >= 9.9_kdp .AND. logp < 10.0_kdp) THEN
     klo = INT(logp*200._kdp) - 1900 ! ... for 8,000 < P < 10,000 bars - step by 0.01 log bars
  ELSEIF (logp >= 10.0_kdp) THEN
     klo = mdimt0 - 1
  END IF
  ! ... adjust index to point to largest pt0 <= px
  ! ...   pt0 - Saturation pressures at selected points
  DO WHILE (px > pt0(klo+1) .AND. klo < mdimt0-1)
     klo = klo + 1
  END DO
  DO WHILE (px < pt0(klo) .AND. klo > 1)
     klo = klo - 1
  END DO
  ! ... Interpolate with the cubic spline
  khi = klo + 1
  h = pt0(khi) - pt0(klo)
  a = (pt0(khi) - px)/h
  b = (px - pt0(klo))/h
  a3 = a**3 - a
  b3 = b**3 - b
  h2 = h**2/6._kdp
  hx = a*t0h(klo) + b*t0h(khi) + (a3*d2t0h(klo) + b3*d2t0h(khi))*h2
  ! ...   t0h - sat'd steam enthalpy at selected points
  ! ...   d2t0h - Sat'd steam enthalpy second derivative w.r.t. pressure
  t0d1 = a*t0d(klo) + b*t0d(khi) + (a3*d2t0d(klo) + b3*d2t0d(khi))*h2
  pt0dph1 = a*pt0dph(klo) + b*pt0dph(khi) +  & 
       (a3*d2pt0dph(klo) + b3*d2pt0dph(khi))*h2
  pt0dpp1 = a*pt0dpp(klo) + b*pt0dpp(khi) +  &
       (a3*d2pt0dpp(klo) + b3*d2pt0dpp(khi))*h2
  t0t1 = a*t0t(klo) + b*t0t(khi) + (a3*d2t0t(klo)+b3*d2t0t(khi))*h2
  pt0tph1 = a*pt0tph(klo) + b*pt0tph(khi) +  &
       (a3*d2pt0tph(klo) + b3*d2pt0tph(khi))*h2
  pt0tpp1 = a*pt0tpp(klo) + b*pt0tpp(khi) +  &
       (a3*d2pt0tpp(klo) + b3*d2pt0tpp(khi))*h2
  t0v1 = a*t0v(klo) + b*t0v(khi) + (a3*d2t0v(klo)+b3*d2t0v(khi))*h2
  pt0vph1 = a*pt0vph(klo) + b*pt0vph(khi) +  &
       (a3*d2pt0vph(klo) + b3*d2pt0vph(khi))*h2
  pt0vpp1 = a*pt0vpp(klo) + b*pt0vpp(khi) +  &
       (a3*d2pt0vpp(klo) + b3*d2pt0vpp(khi))*h2
END SUBROUTINE bndyt0
