SUBROUTINE phsbndy3(px, hw, hs, dhwp, dhsp, swdx, dswdxdp, swvx,  &
     dswvxdp, ssdx, dssdxdp, ssvx, dssvxdp, sstx,  &
     dsstxdp, d2sstxd2p)
  ! ... Purpose:  To calculate the enthaply of sat'd water and sat'd steam
  ! ...         for a given pressure.
  ! ... Note: phsbndy3 is identical to phsbndy1 and phsbndy2 except that
  ! ...         the derivatives w.r.t. P are calculated along the sat'n curves
  ! ...         (i.e. H = HW or HS where HW= f(P).  Also, there is no need
  ! ...         to calculate the derivatives w.r.t. H.
  ! ...         (In phsbndy2, the derivatives w.r.t. P & H calculated
  ! ...            away from the sat'n curves
  ! ... Method: This subroutine uses a second derivative cubic spine.
  ! ... Reference:  Press et al., 1986, Numerical Recipes, the Art of
  ! ...               Scientific Computing, Chapter 3.3 on Cubic Spline
  ! ...               Interpolation, p. 86-87.
  ! ... Range: This spline covers the range in pressure from 0.5 bars
  ! ...           to the critical point, 220.55 bars with 134 selected
  ! ...           pressure-enthalpy pairs.  The pressure enthalpy pairs
  ! ...           were calculated from the steam tables (Haar et al., 1984)
  USE machine_constants, ONLY: kdp
  USE tables
  IMPLICIT NONE
  REAL(kind=kdp), INTENT(IN) :: px
  REAL(kind=kdp), INTENT(OUT) :: hw, hs, dhwp, dhsp, swdx, dswdxdp, swvx, dswvxdp, &
       ssdx, dssdxdp, ssvx, dssvxdp, sstx, dsstxdp, d2sstxd2p
  ! ...   PX - pressure [dyne/cm^2]
  ! ...   HW - enthalpy of sat'd water [erg/g]
  ! ...   HS - enthalpy of sat'd steam [erg/g]
  ! ...   DHWP -  deriv of enthalpy of sat'd water w.r.t P
  ! ...   DHSP -  deriv of enthalpy of sat'd steam w.r.t P
  ! ...   SWDX -  density of sat'd water at P=PX
  ! ...   DSWDXDP -  partial derivative of density of sat'd water at P=PX w.r.t. P
  ! ...   SWVX -  viscosity of sat'd water at P=PX
  ! ...   DSWVXDP -  partial derivative viscosity of sat'd water at P=PX w.r.t. P
  ! ...   SSDX -  density of sat'd steam at P=PX
  ! ...   DSSDXDP -  partial derivative density of sat'd steam at P=PX w.r.t. P
  ! ...   SSVX -  viscsty of sat'd steam at P=PX
  ! ...   DSSVXDP -  partial derivative viscosity of sat'd steam at P=PX w.r.t. P
  ! ...   SSTX -  temp of sat'd steam/water at P=PX
  ! ...   DSSTXDP -  partial derivative temp of sat'd steam/water at P=PX w.r.t. P
  ! ...   D2SSTXD2P -  2nd partial derivative temp of sat'd steam/water at PX wrt P
  !
  INTEGER :: khi, klo
  REAL(kind=kdp) :: a, a2, a3, b, b2, b3, h, h2, pxx
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  !...
  ! ... Special patch for low pressures
  pxx=MAX(px,5.e5_kdp)
  CALL phsbnklo(3,pxx,klo)
  ! ...       check to sure that klo pointer is positioned correctly
  DO WHILE (pxx > ps(klo+1) .AND. klo < mdimws-1)
     klo = klo + 1
  END DO
  DO WHILE (pxx < ps(klo) .AND. klo > 1)
     klo = klo - 1
  END DO
  ! ...     Interpolate with the cubic spline
  khi = klo + 1
  h = ps(khi) - ps(klo)
  a = (ps(khi)-pxx)/h
  b = (pxx-ps(klo))/h
  a3 = a**3 - a
  b3 = b**3 - b
  h2 = h**2/6.0_kdp
  a2 = (-3.0_kdp*a*a+1.0_kdp)*h/6.0_kdp
  b2 = (3.0_kdp*b*b-1.0_kdp)*h/6.0_kdp
  ! ...             water
  hw = a*swh(klo) + b*swh(khi) + (a3*d2swh(klo)+b3*d2swh(khi))*h2
  dhwp = (swh(khi)-swh(klo))/h + a2*d2swh(klo) + b2*d2swh(khi)
  swdx = a*swd(klo) + b*swd(khi) + (a3*d2swd(klo)+b3*d2swd(khi))*h2
  dswdxdp = (swd(khi)-swd(klo))/h + a2*d2swd(klo) + b2*d2swd(khi)
  swvx = a*swv(klo) + b*swv(khi) + (a3*d2swv(klo)+b3*d2swv(khi))*h2
  dswvxdp = (swv(khi)-swv(klo))/h + a2*d2swv(klo) + b2*d2swv(khi)
  ! ...             steam
  hs = a*ssh(klo) + b*ssh(khi) + (a3*d2ssh(klo)+b3*d2ssh(khi))*h2
  dhsp = (ssh(khi)-ssh(klo))/h + a2*d2ssh(klo) + b2*d2ssh(khi)
  ssdx = a*ssd(klo) + b*ssd(khi) + (a3*d2ssd(klo)+b3*d2ssd(khi))*h2
  dssdxdp = (ssd(khi)-ssd(klo))/h + a2*d2ssd(klo) + b2*d2ssd(khi)
  ssvx = a*ssv(klo) + b*ssv(khi) + (a3*d2ssv(klo)+b3*d2ssv(khi))*h2
  dssvxdp = (ssv(khi)-ssv(klo))/h + a2*d2ssv(klo) + b2*d2ssv(khi)
  ! ...          temperature of steam and water
  ! ...          calculated using steam data
  sstx = a*sst(klo) + b*sst(khi) + (a3*d2sst(klo)+b3*d2sst(khi))*h2
  dsstxdp = (sst(khi)-sst(klo))/h + a2*d2sst(klo) + b2*d2sst(khi)
  d2sstxd2p = (ps(khi)-pxx)*d2sst(klo)/h + (pxx-ps(klo))*d2sst(khi)/h
END SUBROUTINE phsbndy3
