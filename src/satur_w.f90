SUBROUTINE satur_w(krtype,pw,hw,sat_w,pswp,pswh)
  ! ... Purpose:  To compute the saturation of water as a function of
  ! ...           capillary pressure and the derivatives w.r.t. p and h
  ! ... Uses the Corey equation that goes with the Corey relative permeability
  ! ...      equation
  ! ... Alternates are a linear equation, a Cooley equation, a van Genuchten equation
  USE machine_constants, ONLY: kdp
  USE parameters, ONLY: patm, pwb, pwr, swr, ssr, bb, cc, gamma, lambda
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: krtype
  REAL(kind=kdp), INTENT(IN) :: pw, hw
  REAL(kind=kdp), INTENT(OUT) :: sat_w, pswp, pswh
  ! ...   KRTYPE (input) - type of pressure-saturation function used
  ! ...            =0 for linear, =1 for Corey
  REAL(kind=kdp) :: pc, pcb, pcr, se, f1, phivg
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.6 $//$Date: 2007/04/10 23:14:34 $'
  !     ------------------------------------------------------------------
  !...
  pc = patm - pw
  pcb = patm - pwb
  IF (krtype == 0) THEN          ! ... linear function
     pcr = patm - pwr
     IF(pc < pcb) THEN
        sat_w = 1._kdp
        pswp = 0._kdp
        pswh = 0._kdp
     ELSE IF(pc > pcr) THEN
        sat_w = swr
        pswp = 0._kdp
        pswh = 0._kdp
     ELSE
        se = 1._kdp - (pc-pcb)/(pcr-pcb)
        sat_w = se*(1._kdp - swr) + swr
        ! ... the derivatives
        pswp = (1._kdp - swr)/(pcr - pcb)
        pswh = 0._kdp        !***...  not sure what this should be
     END IF
  ELSE IF (krtype == 1) THEN      ! ... Corey function
     IF(pc > pcb) THEN
        se = (pcb/pc)**2
        sat_w = se*(1._kdp - swr) + swr
        ! ... the derivatives
        pswp = (1._kdp - swr)*2._kdp*pcb*pcb/(pc*pc*pc)
        pswh = 0._kdp       !***...  not sure what this should be
     ELSE
        sat_w = 1._kdp
        pswp = 0._kdp
        pswh = 0._kdp
     END IF
  ELSEIF(krtype == 3) THEN       ! ... Cooley's modification of Brutsaert function
     IF(pc > pcb) THEN
        se = bb/((pc - pcb)**cc + bb)
        sat_w = se*(1._kdp - swr) + swr
        ! ... the derivatives
        pswp = (1._kdp - swr)*bb*cc*((pc - pcb)**(cc-1._kdp))/ &
             (((pc - pcb)**cc + bb)**2)
        pswh = 0._kdp        !...  not sure what this should be
     ELSE
        sat_w = 1._kdp
        pswp = 0._kdp
        pswh = 0._kdp
     ENDIF
!!$  ELSEIF(krtype == 4) THEN     ! ... van Genuchten function
!!$     IF(pc > 0._kdp) THEN
!!$        phivg = pc/pcb
!!$        f1 = (1._kdp+phivg**lambda) 
!!$        se = f1**(-gamma)
!!$        sat_w = se*(1._kdp - swr) + swr
!!$        ! ... the derivatives
!!$        pswp = -(1._kdp - swr)*((lambda - 1._kdp)*phivg**(lambda-1._kdp))/  &
!!$                  ((f1**(gamma+1._kdp))*pcb)
!!$        pswh = 0._kdp        !...  not sure what this should be
!!$     ELSE
!!$        sat_w = 1._kdp
!!$        pswp = 0._kdp
!!$        pswh = 0._kdp
!!$     ENDIF
  END IF
END SUBROUTINE satur_w
