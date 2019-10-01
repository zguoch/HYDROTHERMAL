SUBROUTINE print_control(privar,time,itime,timchg,timprvar,prvar)
  ! ... Sets the print control flag and stop sign for a given output file
  USE machine_constants, ONLY: kdp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: itime
  REAL(KIND=kdp), INTENT(IN) :: privar, timchg, time
  REAL(KIND=kdp), INTENT(INOUT) :: timprvar
  LOGICAL, INTENT(OUT) :: prvar
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.7 $//$Date: 2003/01/22 22:00:03 $'
  !     ------------------------------------------------------------------
  prvar = .FALSE.
  IF(time >= timchg) THEN
     prvar = .TRUE.
     IF(privar > 0._kdp) timprvar=(1._kdp+INT(time/privar))*privar
  ELSE IF(itime == 0) THEN
     prvar = .TRUE.
  ELSE IF(privar > 0._kdp) THEN
     IF(ABS(timprvar-time) <= 1.e-6_kdp) THEN
        prvar = .TRUE.
        timprvar=(1._kdp+INT(time/privar))*privar
     END IF
  ELSE IF(privar < 0._kdp) THEN
     IF(MOD(itime,INT(ABS(privar))) == 0) prvar = .TRUE.
  END IF
END SUBROUTINE print_control
