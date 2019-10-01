SUBROUTINE table4(hx, px, d0, pd0ph, pd0pp, t0, pt0ph, pt0pp, v0, pv0ph, pv0pp)
  ! ... Purpose:  Lookup table for water properties in the range:
  ! ...         0.5e9 <= H <= 52.0e9 erg/g
  ! ...         240 <= P < 10,000 bar
  ! ... Reference: two-dimensional bicubic interpolation routine is from:
  ! ...      Press et al., 1987, Numerical Recipies, p. 98-100
  USE machine_constants, ONLY: kdp
  USE tables
  IMPLICIT NONE
  ! ...   HX - enthalpy in ergs/g (1 J = 1.0e7 erg)
  ! ...   PX - pressure in dynes/cm^2 (1 bar = 1.0e6 dyne/cm^2)
  ! ...   D0 - density of water g/cm^3
  ! ...      if D0<0 then set D0 to -1*HERGS1(ilo)
  ! ...   PD0PH -  derivative of density w.r.t. enthalpy
  ! ...   PD0PP -  derivative of density w.r.t. pressure
  ! ...   T0 -  temperature of water deg.C
  ! ...   PT0PH -  derivative of temperature w.r.t. enthalpy
  ! ...   PT0PP -  derivative of temperature w.r.t. pressure
  ! ...   V0 -  viscosity of water g/cm-s
  ! ...   PV0PH - derivative of viscosity w.r.t. enthalpy
  ! ...   PV0PP - derivative of viscosity w.r.t. pressure
  REAL(KIND=kdp), INTENT(IN) :: hx, px
  REAL(KIND=kdp), INTENT(OUT) :: d0, pd0ph, pd0pp, t0, pt0ph, pt0pp, v0, pv0ph, pv0pp
  !
  REAL(KIND=kdp) :: ansy, ansy1, ansy2,  tt, uu
  INTEGER :: i, ihi, ilo, jhi, jlo, l
  REAL(KIND=kdp), DIMENSION(3) :: u, puph, pupp
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3 $//$Date: 2005/12/02 23:35:20 $'
  ! ... ------------------------------------------------------------------
  ! ... position pointer for enthalpy - ilo
  ilo = 1
  IF (hx >= 52.0D9) THEN
     ilo = idim4 - 1
     ! ... for 40.0 <= H < 52.0e9, H steps by 4.0e9 ergs/g
  ELSE IF (hx >= 40.0e9_kdp) THEN
     ilo = IDINT(hx*1.0e-9_kdp/4._kdp) + 25
     ! ... for 10.0 <= H < 40.0e9, H steps by 2.0e9 ergs/g
  ELSE IF (hx >= 10.0e9_kdp) THEN
     ilo = IDINT(hx*1.0e-9_kdp/2._kdp) + 15
     ! ... for H=0.5e9 to 10.0e9 ergs/g step by 0.5e9 ergs/g
  ELSE IF (hx >= 0.5e9_kdp) THEN
     ilo = IDINT(hx*1.0e-8_kdp/5._kdp)
  END IF
  ! ... check position of pointer - ilo
  DO WHILE (hx > hergs4(ilo+1) .AND. ilo < idim4-1)
     ilo = ilo + 1
  END DO
  DO WHILE (hx < hergs4(ilo) .AND. ilo > 1)
     ilo = ilo - 1
  END DO
  ihi = ilo + 1
  ! ... position pointer for pressure - jlo
  ! ... if 240 <= P < 10,000 bars step by 0.05 log bars
  ! ...     IF (PX.GE.240.0D6.AND.PX.LT.1.0D10)
  jlo = IDINT(LOG10(px)*20._kdp) - 166
  IF (px >= 1.0e10_kdp) jlo = jdim4 - 1
  ! ... check position of pointer - jlo
  DO WHILE (px > pdyne4(jlo+1) .AND. jlo < jdim4-1)
     jlo = jlo + 1
  END DO
  DO WHILE (px < pdyne4(jlo) .AND. jlo > 1)
     jlo = jlo - 1
  END DO
  jhi = jlo + 1
  tt = (hx-hergs4(ilo))/(hergs4(ihi)-hergs4(ilo))
  uu = (px-pdyne4(jlo))/(pdyne4(jhi)-pdyne4(jlo))
  ! ... L=1:Density, L=2:Temperature, L=3:Viscosity
  DO  l=1,3
     ansy = 0._kdp
     ansy2 = 0._kdp
     ansy1 = 0._kdp
     DO  i=4,1,-1
        ansy = tt*ansy + ((ctbl4(ilo,jlo,l,i,4)*uu+ctbl4(ilo,jlo,l,i,3))  &
             *uu+ctbl4(ilo,jlo,l,i,2))*uu + ctbl4(ilo,jlo,l,i,1)
        ansy2 = tt*ansy2 +  &
             (3._kdp*ctbl4(ilo,jlo,l,i,4)*uu+2._kdp*ctbl4(ilo,jlo,l,  &
             i,3))*uu + ctbl4(ilo,jlo,l,i,2)
        ansy1 = uu*ansy1 +  &
             (3._kdp*ctbl4(ilo,jlo,l,4,i)*tt+2._kdp*ctbl4(ilo,jlo,l,  &
             3,i))*tt + ctbl4(ilo,jlo,l,2,i)
     END DO
     u(l) = ansy
     puph(l) = ansy1/(hergs4(ihi)-hergs4(ilo))
     pupp(l) = ansy2/(pdyne4(jhi)-pdyne4(jlo))
  END DO
  d0 = u(1)
  ! ... if point HX,PX is between the tabulated values and the T=0 curve
  ! ...    set the value of D0 to -1* enthalpy of the next (higher)
  ! ...    tabulated value for H.
  IF (d0 <= 0._kdp) d0 = -1._kdp*(hergs4(ihi)+0.0001e9_kdp)
  pd0ph = puph(1)
  pd0pp = pupp(1)
  t0 = u(2)
  pt0ph = puph(2)
  pt0pp = pupp(2)
  v0 = u(3)
  pv0ph = puph(3)
  pv0pp = pupp(3)
END SUBROUTINE table4
