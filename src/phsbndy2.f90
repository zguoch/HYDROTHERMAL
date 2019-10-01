SUBROUTINE phsbndy2(px, hw, hs, dhwp, dhsp, swdx, pswdxph,  &
     pswdxpp, swtx, pswtxph, pswtxpp, swvx,  &
     pswvxph, pswvxpp, ssdx, pssdxph, pssdxpp,  &
     sstx, psstxph, psstxpp, ssvx, pssvxph,  &
     pssvxpp)
  !     Purpose:  To calculate the enthaply of sat'd water and sat'd steam
  !             for a given pressure.
  !     Note: phsbndy2 is identical to phsbndy1 except that fewer parameters
  !             are calculated in phsbndy1 (only HW and HS).
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
  REAL(kind=kdp), INTENT(OUT) :: hw, hs, dhwp, dhsp, swdx, pswdxph, pswdxpp, swtx, &
       pswtxph, pswtxpp, swvx, pswvxph, pswvxpp, ssdx, pssdxph, pssdxpp, sstx, psstxph, &
       psstxpp, ssvx, pssvxph, pssvxpp
  !       PX - pressure in dynes/cm^2
  !       HW - enthalpy of sat'd water in erg/g
  !       HS - enthalpy of sat'd steam in erg/g
  !       DHWP -  deriv of enthalpy of sat'd water w.r.t P
  !       DHSP -  deriv of enthalpy of sat'd steam w.r.t P
  !       SWDX -  density of sat'd water at P=PX
  !       PSWDXPH -  partial deriv density of sat'd water at P=PX w.r.t. H
  !       PSWDXPP -  partial deriv density of sat'd water at P=PX w.r.t. P
  !       SWTX -  temp of sat'd water at P=PX
  !       PSWTXPH -  partial deriv temp of sat'd water at P=PX w.r.t. H
  !       PSWTXPP -  partial deriv temp of sat'd water at P=PX w.r.t. P
  !       SWVX -  viscsty of sat'd water at P=PX
  !       PSWVXPH -  partial deriv viscsty of sat'd water at P=PX w.r.t. H
  !       PSWVXPP -  partial deriv viscsty of sat'd water at P=PX w.r.t. P
  !       SSDX -  density of sat'd steam at P=PX
  !       PSSDXPH -  partial deriv density of sat'd steam at P=PX w.r.t. H
  !       PSSDXPP -  partial deriv density of sat'd steam at P=PX w.r.t. P
  !       SSTX -  temp of sat'd steam at P=PX
  !       PSSTXPH -  partial deriv temp of sat'd steam at P=PX w.r.t. H
  !       PSSTXPP -  partial deriv temp of sat'd steam at P=PX w.r.t. P
  !       SSVX -  viscsty of sat'd steam at P=PX
  !       PSSVXPH -  partial deriv viscsty of sat'd steam at P=PX w.r.t. H
  !       PSSVXPP -  partial deriv viscsty of sat'd steam at P=PX w.r.t. P
  !
  INTEGER :: khi, klo
  REAL(kind=kdp) :: a, a3, b, b3, h, h2
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.2 $//$Date: 2000/12/12 21:27:01 $'
  !     ------------------------------------------------------------------
  !...
  CALL phsbnklo (2, px, klo)
  !           check to sure that klo pointer is positioned correctly
  DO WHILE (px > ps(klo+1) .AND. klo < mdimws-1)
     klo = klo + 1
  END DO
  DO WHILE (px < ps(klo) .AND. klo > 1)
     klo = klo - 1
  END DO
  !         Interpolate with the cubic spline
  khi = klo + 1
  h = ps(khi) - ps(klo)
  a = (ps(khi)-px)/h
  b = (px-ps(klo))/h
  a3 = a**3 - a
  b3 = b**3 - b
  h2 = h**2/6.0_kdp
  !                 water
  hw = a*swh(klo) + b*swh(khi) + (a3*d2swh(klo)+b3*d2swh(khi))*h2
  dhwp = (swh(khi)-swh(klo))/h + (-3.0_KDP*a*a+1.0_KDP)*d2swh(klo)*h/6.0_KDP  &
       + (3.0_KDP*b*b-1.0_KDP)*d2swh(khi)*h/6.0_KDP
  swdx = a*swd(klo) + b*swd(khi) + (a3*d2swd(klo)+b3*d2swd(khi))*h2
  pswdxph = a*pswdph(klo) + b*pswdph(khi)  &
       + (a3*d2pswdph(klo)+b3*d2pswdph(khi))*h2
  pswdxpp = a*pswdpp(klo) + b*pswdpp(khi)  &
       + (a3*d2pswdpp(klo)+b3*d2pswdpp(khi))*h2
  swtx = a*swt(klo) + b*swt(khi) + (a3*d2swt(klo)+b3*d2swt(khi))*h2
  pswtxph = a*pswtph(klo) + b*pswtph(khi)  &
       + (a3*d2pswtph(klo)+b3*d2pswtph(khi))*h2
  pswtxpp = a*pswtpp(klo) + b*pswtpp(khi)  &
       + (a3*d2pswtpp(klo)+b3*d2pswtpp(khi))*h2
  swvx = a*swv(klo) + b*swv(khi) + (a3*d2swv(klo)+b3*d2swv(khi))*h2
  pswvxph = a*pswvph(klo) + b*pswvph(khi)  &
       + (a3*d2pswvph(klo)+b3*d2pswvph(khi))*h2
  pswvxpp = a*pswvpp(klo) + b*pswvpp(khi)  &
       + (a3*d2pswvpp(klo)+b3*d2pswvpp(khi))*h2
  !                 steam
  hs = a*ssh(klo) + b*ssh(khi) + (a3*d2ssh(klo)+b3*d2ssh(khi))*h2
  dhsp = (ssh(khi)-ssh(klo))/h + (-3.0_KDP*a*a+1.0_KDP)*d2ssh(klo)*h/6.0_KDP  &
       + (3.0_KDP*b*b-1.0_KDP)*d2ssh(khi)*h/6.0_KDP
  ssdx = a*ssd(klo) + b*ssd(khi) + (a3*d2ssd(klo)+b3*d2ssd(khi))*h2
  pssdxph = a*pssdph(klo) + b*pssdph(khi)  &
       + (a3*d2pssdph(klo)+b3*d2pssdph(khi))*h2
  pssdxpp = a*pssdpp(klo) + b*pssdpp(khi)  &
       + (a3*d2pssdpp(klo)+b3*d2pssdpp(khi))*h2
  sstx = a*sst(klo) + b*sst(khi) + (a3*d2sst(klo)+b3*d2sst(khi))*h2
  psstxph = a*psstph(klo) + b*psstph(khi)  &
       + (a3*d2psstph(klo)+b3*d2psstph(khi))*h2
  psstxpp = a*psstpp(klo) + b*psstpp(khi)  &
       + (a3*d2psstpp(klo)+b3*d2psstpp(khi))*h2
  ssvx = a*ssv(klo) + b*ssv(khi) + (a3*d2ssv(klo)+b3*d2ssv(khi))*h2
  pssvxph = a*pssvph(klo) + b*pssvph(khi)  &
       + (a3*d2pssvph(klo)+b3*d2pssvph(khi))*h2
  pssvxpp = a*pssvpp(klo) + b*pssvpp(khi)  &
       + (a3*d2pssvpp(klo)+b3*d2pssvpp(khi))*h2
END SUBROUTINE phsbndy2
