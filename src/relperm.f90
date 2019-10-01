SUBROUTINE relperm(krtype, pw, sw, pswp, pswh, xkw, pxkwp,  &
     pxkwh, xks, pxksp, pxksh, ind)
  ! ... Purpose:  To calculate the relative permeability of water and steam
  ! ...           using either the liner, corey, or fracture functions if water
  ! ...           sat'n is between residual water sat'n & 1.0 - residual
  ! ...           steam sat'n
  !
  ! ... Reference:  Ingebritsen, S.E., 1986, Vapor-dominated zones within
  ! ...      hydrothermal convection systems: evolution and natural state:
  ! ...      Unpub. Ph.D. thesis, Stanford University, p. 150-152.
  USE machine_constants, ONLY: kdp
  USE f_units
  USE control
  USE parameters, ONLY: patm, pwb, pwr, swr, ssr, bb, cc, gamma, lambda
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: krtype
  REAL(kind=kdp), INTENT(IN) :: pw, sw, pswp, pswh
  REAL(kind=kdp), INTENT(OUT) :: xkw, pxkwp, pxkwh, xks, pxksp, pxksh
  INTEGER, INTENT(IN) :: ind
  ! ... Passed Variables:
  ! ...   KRTYPE (input) - type of relative permeability function used
  ! ...            =0 for linear, =1 for Corey, =2 for fracture
  ! ...            =3 for Cooley, =4 for van Genuchten-Mualem
  ! ...   SW (input) - saturation of water (1.0 = sat'd)
  ! ...   PSWP (input) - derivative of the saturation of water w.r.t. P
  ! ...   PSWH (input) - derivative of the saturation of water w.r.t. H
  ! ...   SWR (input) - residual water saturation
  ! ...   SSR (input) - residual steam saturation
  ! ...   XKW (retn) - relative permeability of water
  ! ...   PXKWP (retn) - deriv. of relative permeability of water w.r.t. P
  ! ...   PXKWH (retn) - deriv. of relative permeability of water w.r.t. H
  ! ...   XKS (retn) - relative permeability of steam
  ! ...   PXKSP (retn) - deriv. of relative permeability of steam w.r.t. P
  ! ...   PXKSH (retn) - deriv. of relative permeability of steam w.r.t. H
  !
  REAL(kind=kdp) :: sss, sssu, lm1, f1, f2, phivg, pc, pcb
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5 $//$Date: 2006/10/03 23:29:10 $'
  !     ------------------------------------------------------------------
  ! ... Check if water saturation is greater than 1.0 - residual steam saturation
  IF (sw >= (1._kdp-ssr)) THEN
     !*****special patch for fast draining to steady state
!!$***  IF (sw >= (1._kdp-ssr) .or. ind == 5) THEN
     xkw = 1._kdp
     xks = 0._kdp
     pxkwp = 0._kdp
     pxkwh = 0._kdp
     pxksp = 0._kdp
     pxksh = 0._kdp
     ! ... Check if water saturation is less than residual water saturation
  ELSE IF (sw <= swr) THEN
     xkw = 0._kdp
     xks = 1._kdp
     pxkwp = 0._kdp
     pxkwh = 0._kdp
     pxksp = 0._kdp
     pxksh = 0._kdp
  ELSE IF (krtype == 0) THEN
     ! ... linear
     xkw = (sw-swr)/((1._kdp-swr)-ssr)
     xkw = MIN(1._kdp, xkw)
     xks = 1._kdp - xkw
     pxkwp = pswp/((1._kdp-swr)-ssr)
     pxkwh = pswh/((1._kdp-swr)-ssr)
     pxksp = -pxkwp
     pxksh = -pxkwh
  ELSE IF (krtype == 1.OR. krtype == 3) THEN
     ! ... Corey also used with Cooley Sw(p) curve
     sssu = (1._kdp-swr) - ssr
     sss = (sw-swr)/sssu
     xkw = sss**4
     pxkwp = 4._kdp*(sss**3)*pswp/sssu
     pxkwh = 4._kdp*(sss**3)*pswh/sssu
     IF (xkw >= 1._kdp) THEN
        xkw = 1._kdp
        xks = 0._kdp
     ELSE
        xks = (1._kdp-sss**2)*(1._kdp-sss)**2
        xks = MIN(1._kdp, xks)
        pxksp = -2._kdp*sss*pswp/sssu*((1._kdp-sss)**2)  &
             + (1._kdp-sss**2)*2._kdp*(1._kdp-sss)*(-pswp)/sssu
        pxksh = -2._kdp*sss*pswh/sssu*((1._kdp-sss)**2)  &
             + (1._kdp-sss**2)*2._kdp*(1._kdp-sss)*(-pswh)/sssu
     END IF
  ELSE IF (krtype == 2) THEN
     ! ... fracture flow
     sssu = (1._kdp-swr) - ssr
     sss = (sw-swr)/sssu
     xkw = sss**4
     xkw = MIN(1._kdp, xkw)
     xks = 1._kdp - xkw
     pxkwp = 4._kdp*(sss**3)*pswp/sssu
     pxkwh = 4._kdp*(sss**3)*pswh/sssu
     pxksp = -pxkwp
     pxksh = -pxkwh
!!$  ELSE IF (krtype == 4) THEN
!!$     !---van Genuchten-Mualem
!!$     pc = patm - pw
!!$     pcb = patm - pwb
!!$     sssu = (1._kdp-swr) - ssr
!!$     sss = (sw-swr)/sssu
!!$     xkw = sqrt(sss)*(1._kdp - (1._kdp - sss**(1._kdp/gamma))**gamma)**2
!!$     xkw = MIN(1._kdp, xkw)
!!$     xks = 1._kdp - xkw
!!$     lm1 = lambda-1._kdp
!!$     phivg = pc/pcb
!!$     f1 = (1._kdp+phivg**lambda) 
!!$     pxkwp = lm1*(phivg**lm1)*(f1**(-2.5_kdp*gamma))*(f1**gamma - phivg**(lambda*gamma))*  &
!!$          (phivg**(lambda*(gamma-1._kdp))*((5._kdp*(phivg**lambda))/(2._kdp*f1) - 2._kdp) -  &
!!$          0.5_kdp*f1**(gamma-1._kdp))*(1._kdp/pcb)
!!$     f2 = (1._kdp - ( 1._kdp - sss**(1._kdp/gamma))**gamma)
!!$     pxkwh = sqrt(sss)*f2*(0.5_kdp*f2 - 2._kdp*(1._kdp - sss**(1._kdp/gamma))**(gamma-1._kdp))*  &
!!$          (1._kdp/(1._kdp - swr))*pswh
!!$     pxksp = -pxkwp
!!$     pxksh = -pxkwh
  ELSE
     ! ...  KRTYPE is not 0,1,2,3 or 4--stop
     WRITE (fustdout, 9005) krtype
     WRITE (fuclog, 9005) krtype
9005 FORMAT (' *** STOP ***   ERROR with relative permeability type: ',i2)
     ierr(51) = .TRUE.
  END IF
END SUBROUTINE relperm
