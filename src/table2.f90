SUBROUTINE table2(hx, px, d0, pd0ph, pd0pp, t0, pt0ph, pt0pp, v0, pv0ph, pv0pp)
  ! ... Purpose:  Lookup table for water properties in the range:
  ! ...         26.0D9 <= H <= 52.0D9 erg/g
  ! ...         0.5 <= P < 240 bar
  ! ... Reference: two-dimensional bicubic interpolation routine is from
  ! ...      Press et al., 1987, Numerical Recipies, p. 98-100
  USE machine_constants, ONLY: kdp
  USE tables
  IMPLICIT NONE
  ! ...   HX - enthalpy in ergs/g (1 J = 1.0D7 erg)
  ! ...   PX - pressure in dynes/cm^2 (1 bar = 1.0D6 dyne/cm^2)
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
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4 $//$Date: 2006/09/06 21:04:28 $'
  ! ... ------------------------------------------------------------------
  ! ... position pointer for enthalpy - ilo
  ilo = 1
  IF (hx >= 26.0e9_kdp .AND. hx < 30.0e9_kdp) THEN  ! ... 26.e9 < H < 30.e9, H steps by 0.5e9 ergs/g
     ilo = IDINT(hx*1.0e-8_kdp/5._kdp) - 51
  ELSE IF (hx >= 30.0e9_kdp .AND. hx < 40.0e9_kdp) THEN     ! ... for 30.e9 < H < 40.e9, H steps by 2.e9 ergs/g
     ilo = IDINT(hx*1.0e-9_kdp/2._kdp) - 6
  ELSE IF (hx >= 40.0e9_kdp .AND. hx < 52.0e9_kdp) THEN     ! ... for 40.e9 < H < 52.e9, H steps by 5.e9 ergs/g
     ilo = IDINT(hx*1.0e-9_kdp/4._kdp) + 4
  END IF
  ! ... check position of pointer - ilo
  DO WHILE (hx > hergs2(ilo+1) .AND. ilo < idim2-1)
     ilo = ilo + 1
  END DO
  DO WHILE (hx < hergs2(ilo) .AND. ilo > 1)
     ilo = ilo - 1
  END DO
  ihi = ilo + 1
  ! ... position pointer for pressure - jlo
  ! ... if P < 5 bars
  IF (px < 1.0e6_kdp) THEN
     jlo = 1
  ELSE IF (px < 3.0e6_kdp) THEN
     jlo = 2
  ELSE IF (px < 5.0e6_kdp) THEN
     jlo = 3
  ELSE IF (px < 100.0e6_kdp) THEN         ! ... if 5<P<100 bars step by 5 bars
     jlo = IDINT(px*1.0e-6_kdp/5._kdp) + 3
  ELSE                                    ! ... if 100<P<240 bars step by 10 bars
     jlo = IDINT(px*1.0e-6_kdp/10._kdp) + 13
  END IF
  ! ... check position of pointer - jlo
  DO WHILE (px > pdyne2(jlo+1) .AND. jlo < jdim2-1)
     jlo = jlo + 1
  END DO
  DO WHILE (px < pdyne2(jlo) .AND. jlo > 1)
     jlo = jlo - 1
  END DO
  jhi = jlo + 1
  tt = (hx-hergs2(ilo))/(hergs2(ihi)-hergs2(ilo))
  uu = (px-pdyne2(jlo))/(pdyne2(jhi)-pdyne2(jlo))
  ! ... l=1:Density, l=2:Temperature, l=3:Viscosity
  DO  l=1,3
     ansy = 0._kdp
     ansy2 = 0._kdp
     ansy1 = 0._kdp
     DO  i=4,1,-1
        ansy = tt*ansy + ((ctbl2(ilo,jlo,l,i,4)*uu+ctbl2(ilo,jlo,l,i,3))  &
             *uu+ctbl2(ilo,jlo,l,i,2))*uu + ctbl2(ilo,jlo,l,i,1)
        ansy2 = tt*ansy2 +  &
             (3._kdp*ctbl2(ilo,jlo,l,i,4)*uu+2._kdp*ctbl2(ilo,jlo,l,i,3))  &
             *uu + ctbl2(ilo,jlo,l,i,2)
        ansy1 = uu*ansy1 +  &
             (3._kdp*ctbl2(ilo,jlo,l,4,i)*tt+2._kdp*ctbl2(ilo,jlo,l,3,i))  &
             *tt + ctbl2(ilo,jlo,l,2,i)
     END DO
     u(l) = ansy
     puph(l) = ansy1/(hergs2(ihi)-hergs2(ilo))
     pupp(l) = ansy2/(pdyne2(jhi)-pdyne2(jlo))
  END DO
  d0 = u(1)
  ! ...    if point HX,PX is between the tabulated values and the two-phase curve
  ! ...      set the value of D0 to -1* enthalpy of the next (higher)
  ! ...      tabulated value
  IF (d0 <= 0._kdp) d0 = -1._kdp*(hergs2(ihi)+0.0001e9_kdp)
  pd0ph = puph(1)
  pd0pp = pupp(1)
  t0 = u(2)
  pt0ph = puph(2)
  pt0pp = pupp(2)
  v0 = u(3)
  pv0ph = puph(3)
  pv0pp = pupp(3)
END SUBROUTINE table2
