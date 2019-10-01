SUBROUTINE phsbnklo(isub, px, klo)
  !     Purpose:  To estimate position in the lookup table for a given pressure.
  !     Range: This spline covers the range in pressure from 0.5 bar
  !               to the critical point, 220.55 bar, with 134 
  !               pressure-enthalpy data pairs.  The pressure-enthalpy pairs
  !               were calculated from the steam tables (Haar et al., 1984)
  USE machine_constants, ONLY: kdp
  USE f_units
  USE control
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: isub
  REAL(kind=kdp), INTENT(IN):: px
  INTEGER, INTENT(OUT) :: klo
  !       ISUB - number of phsbndy routine that made the call
  !       PX   - pressure [dynes/cm^2]
  !       KLO  - index in lookup table
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  !...
  IF (px < 0.5e6_kdp .OR. px >= 220.55e6_kdp) THEN
     WRITE(fustdout,9005) isub, px
     WRITE(fuclog,9005) isub, px
9005 FORMAT (  &
          ' ***ERROR:Pressure beyond limits of cubic spline in '//  &
          ' phsbnklo'/tr10, 'Called by PHSBNDY',i1,' P=',1PG14.6)
     ierr(72) = .TRUE.
     !..        return
  END IF
  ! ... estimate position in the lookup table
  IF (px < 210.0e6_kdp) THEN
     ! ... if pressure is < 0.6 bar
     IF (px < 0.6e6_kdp) THEN
        klo = 1
        ! ... for 0.6 bar < PX < 3 bar   PS steps by 0.1 bar
     ELSE IF (px < 3.0e6_kdp) THEN
        klo = IDINT(px*1.0e-5_kdp) - 4
        ! ... for 3 bar < PX < 30 bar    PS steps by 1 bar
     ELSE IF (px < 30.0E6_KDP) THEN
        klo = IDINT(px*1.0e-6_kdp) + 23
        ! ... for 30 bar < PX < 210 bar  PS steps by 5 bar
     ELSE
        klo = IDINT(DBLE(IDINT(px*1.0e-6_kdp))/5._kdp) + 47
     END IF
     ! ... for 210 bar < PX < 215 bar    PS steps by 1 bar
  ELSE IF (px < 215.0e6_kdp) THEN
     klo = IDINT(px*1.0e-6_kdp) - 121
     ! ... for 215 bar < PX < 220 bar    PS steps by 0.2 bar
  ELSE IF (px < 220.0e6_kdp) THEN
     klo = IDINT(DBLE(IDINT(px*1.0e-5_kdp)-2150)/2.0_kdp) + 94
     ! ... for 220 bar < PX < 220.50 bar PS steps by 0.05 bar
  ELSE IF (px < 220.50E6_KDP) THEN
     klo = IDINT(DBLE(IDINT(px*1.0e-4_kdp)-22000)/5._kdp) + 119
     ! ... for 220.50 bar < PX < 220.55 bar PS steps by 0.01 bar
  ELSE
     klo = IDINT(px*1.0e-4_kdp) - 22050 + 129
  END IF
END SUBROUTINE phsbnklo
