SUBROUTINE timestep_criterion(ltmcntl,atmcntl)
  !     Purpose: print the reason for the time step change
  USE f_units
  USE control, ONLY: itmcntl, ierr, ioptpr
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ltmcntl
  CHARACTER (LEN=1), INTENT(OUT) :: atmcntl
  !
  CHARACTER (LEN=80) :: string
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4 $//$Date: 2004/09/20 21:51:00 $'
  !     ------------------------------------------------------------------
  IF (ltmcntl == 0) THEN     ! ... not adjusted
     string = 'Time step not adjusted' 
     atmcntl = 'X'
  ELSE IF (ltmcntl == 1) THEN     ! ... water saturation change
     itmcntl(1) = itmcntl(1) + 1
     atmcntl = 'W'
     string = 'Time step set for maximum Water Saturation change'
  ELSE IF (ltmcntl == 2) THEN     ! ... maximum pressure change
     itmcntl(2) = itmcntl(2) + 1
     atmcntl = 'P'
     string = 'Time step set for maximum Pressure change'
  ELSE IF (ltmcntl == 3) THEN     ! ... maximum enthalpy change
     itmcntl(3) = itmcntl(3) + 1
     atmcntl = 'H'
     string = 'Time step set for maximum Enthalpy change'
  ELSE IF (ltmcntl == 4) THEN     ! ... maximum time step increase
     itmcntl(4) = itmcntl(4) + 1
     atmcntl = 'M'
     string = 'Time step set by factor for maximum time step increase'
  ELSE IF (ltmcntl == 5) THEN     ! ... reset to value before writing
     itmcntl(5) = itmcntl(5) + 1
     atmcntl = 'R'
     string = 'Time step reset to value before writing output'
  ELSE IF (ltmcntl == 6) THEN     ! ... initial time step length
     itmcntl(6) = itmcntl(6) + 1
     atmcntl = 'I'
     string = 'Time step set to input value for initial step'
  ELSE IF (ltmcntl == 7) THEN     ! ... read new input data
     itmcntl(7) = itmcntl(7) + 1
     atmcntl = 'D'
     string = 'Time step set to time for reading new input data'
  ELSE IF (ltmcntl == 8) THEN     ! ... print time
     itmcntl(8) = itmcntl(8) + 1
     atmcntl = 'O'
     string = 'Time step set to time for writing output data'
  ELSE IF (ltmcntl == 9) THEN     ! ... before read time
     itmcntl(9) = itmcntl(9) + 1
     atmcntl = 'R'
     string = 'Time step reset to value before reading new data'
  ELSE IF (ltmcntl == 10) THEN     ! ... maximum time step
     itmcntl(10) = itmcntl(10) + 1
     atmcntl = 'M'
     string = 'Time step set to maximum allowed value'
  ELSE                          ! ... error in value of ltmcntl
     WRITE (string,9040) '**** Invalid time step control flag ****   LTMCNTL:',ltmcntl
9040 FORMAT (a,i3)
     ierr(125) = .TRUE.
     atmcntl = ' '
  END IF
  IF (ioptpr(3) > 0) THEN
     WRITE (fustdout,'(tr8,a)') TRIM(string)
     WRITE (fuclog,'(tr8,a)') TRIM(string)
  END IF
END SUBROUTINE timestep_criterion
