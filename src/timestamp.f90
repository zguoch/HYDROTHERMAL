SUBROUTINE timestamp(date_time)
  CHARACTER(LEN=*) :: date_time
  CHARACTER(LEN=3), DIMENSION(12), PARAMETER :: months=(/'Jan', 'Feb', 'Mar', 'Apr',  &
       'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/)
  INTEGER, DIMENSION(8) :: elements
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5 $//$Date: 2003/10/15 22:03:11 $'
  !     ------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=elements)
  IF(elements(1) /= -HUGE(0)) THEN
!!$     IF(elements(1) < 2000) THEN
!!$        elements(1) = elements(1) - 1900
!!$     ELSE
!!$        elements(1) = elements(1) - 2000
!!$     ENDIF
     WRITE(date_time,1001) elements(3),months(elements(2)),elements(1),  &
          elements(5),elements(6),elements(7)
1001 FORMAT(I2.2,TR1,A3,TR1,I4,TR1,I2.2,':',I2.2,':',I2.2)
  ELSE
     date_time = ' '
  ENDIF
END SUBROUTINE timestamp
