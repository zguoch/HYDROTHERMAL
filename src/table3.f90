SUBROUTINE table3(hx, px, d0, pd0ph, pd0pp, t0, pt0ph, pt0pp, v0, pv0ph, pv0pp)
  ! ... Purpose:  Lookup table for water properties in the range:
  ! ...         16.0e9 <= H <= 26.0e9 ergs/g
  ! ...         150 <= P < 240 bars
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
  ! ... for 16.0 < H < 20.5e9, H steps by 0.5e9 ergs/g
  IF (hx >= 16.0e9_kdp .AND. hx < 20.5e9_kdp) THEN
     ilo = IDINT(hx*1.0e-8_kdp/5._kdp) - 31
     ! ... for 20.5 < H < 21.1e9
  ELSE IF (hx >= 20.5e9_kdp .AND. hx < 21.5e9_kdp) THEN
     IF (hx >= 21.1e9_kdp) THEN
        ilo = 17
     ELSE IF (hx >= 21.0e9_kdp) THEN
        ilo = 16
     ELSE IF (hx >= 20.9e9_kdp) THEN
        ilo = 15
     ELSE IF (hx >= 20.86e9_kdp) THEN
        ilo = 14
     ELSE IF (hx >= 20.8e9_kdp) THEN
        ilo = 13
     ELSE IF (hx >= 20.7e9_kdp) THEN
        ilo = 12
     ELSE IF (hx >= 20.6e9_kdp) THEN
        ilo = 11
     ELSE
        ilo = 10
     END IF
     ! ... for 21.5 < H < 26.0e9, H steps by 0.5e9 ergs/g
  ELSE IF (hx >= 21.5e9_kdp .AND. hx < 26.0e9_kdp) THEN
     ilo = IDINT(hx*1.0e-8_kdp/5._kdp) - 25
  END IF
  ! ... check position of pointer - ilo
  DO WHILE (hx > hergs3(ilo+1) .AND. ilo < idim3-1)
     ilo = ilo + 1
  END DO
  DO WHILE (hx < hergs3(ilo) .AND. ilo > 1)
     ilo = ilo - 1
  END DO
  ihi = ilo + 1
  ! ... position pointer for pressure - jlo
  ! ... if 150 < P< 210 bars step by 5 bars
  IF (px < 150.0e6_kdp) THEN
     jlo = 1
  ELSE IF (px < 210.0e6_kdp) THEN
     jlo = IDINT(px*1.0e-6_kdp/5._kdp) - 29
     ! ... if 210 < P< 220 bars step by 1 bars
  ELSE IF (px < 220.0e6_kdp) THEN
     jlo = IDINT(px*1.0e-6_kdp) - 197
     ! ... if 220 < P< 220.5 bars step by 0.1 bars
  ELSE IF (px < 220.5e6_kdp) THEN
     jlo = IDINT(px*1.0e-5_kdp) - 2177
  ELSE IF (px < 220.55e6_kdp) THEN
     jlo = 28
  ELSE IF (px < 220.6e6_kdp) THEN
     jlo = 29
     ! ... if 220.6 < P< 221.0 bars step by 0.1 bars
  ELSE IF (px < 221.0e6_kdp) THEN
     jlo = IDINT(px*1.0e-5_kdp) - 2176
     ! ... if 221 < P< 230 bars step by 1 bars
  ELSE IF (px < 230.0e6_kdp) THEN
     jlo = IDINT(px*1.0e-6_kdp) - 187
  ELSE IF (px < 235.0e6_kdp) THEN
     jlo = 43
  ELSE
     jlo = 44
  END IF
  ! ... check position of pointer - jlo
  DO WHILE (px > pdyne3(jlo+1) .AND. jlo < jdim3-1)
     jlo = jlo + 1
  END DO
  DO WHILE (px < pdyne3(jlo) .AND. jlo > 1)
     jlo = jlo - 1
  END DO
  jhi = jlo + 1
  tt = (hx-hergs3(ilo))/(hergs3(ihi)-hergs3(ilo))
  uu = (px-pdyne3(jlo))/(pdyne3(jhi)-pdyne3(jlo))
  ! ... L=1:Density, L=2:Temperature, L=3:Viscosity
  DO  l=1,3
     ansy = 0._kdp
     ansy2 = 0._kdp
     ansy1 = 0._kdp
     DO  i=4,1,-1
        ansy = tt*ansy + ((ctbl3(ilo,jlo,l,i,4)*uu+ctbl3(ilo,jlo,l,i,3))  &
             *uu+ctbl3(ilo,jlo,l,i,2))*uu + ctbl3(ilo,jlo,l,i,1)
        ansy2 = tt*ansy2 +  &
             (3._kdp*ctbl3(ilo,jlo,l,i,4)*uu+2._kdp*ctbl3(ilo,jlo,l,  &
             i,3))*uu + ctbl3(ilo,jlo,l,i,2)
        ansy1 = uu*ansy1 +  &
             (3._kdp*ctbl3(ilo,jlo,l,4,i)*tt+2._kdp*ctbl3(ilo,jlo,l,  &
             3,i))*tt + ctbl3(ilo,jlo,l,2,i)
     END DO
     u(l) = ansy
     puph(l) = ansy1/(hergs3(ihi)-hergs3(ilo))
     pupp(l) = ansy2/(pdyne3(jhi)-pdyne3(jlo))
  END DO
  d0 = u(1)
  ! ... if point HX,PX is between the tabulated values and the two-phase curve
  ! ...    set the value of D0 to -1* enthalpy of the next lower
  ! ...    tabulated value for H < H critical, and
  ! ...    set the value of D0 to -1* enthalpy of the next (higher)
  ! ...    tabulated value for H > H critical.
  IF (d0 <= 0._kdp .AND. hx <= 20.86e9_kdp) d0 = -1._kdp*(hergs3(ilo)-0.0001e9_kdp)
  IF (d0 <= 0._kdp .AND. hx > 20.86e9_kdp) d0 = -1._kdp*(hergs3(ihi)+0.0001e9_kdp)
  pd0ph = puph(1)
  pd0pp = pupp(1)
  t0 = u(2)
  pt0ph = puph(2)
  pt0pp = pupp(2)
  v0 = u(3)
  pv0ph = puph(3)
  pv0pp = pupp(3)
END SUBROUTINE table3
