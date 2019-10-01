FUNCTION uppercase(string) RESULT(outstring)
  ! ... Change string to uppercase
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: string
  CHARACTER(LEN=len(string)) :: outstring
  !
  CHARACTER(LEN=1) :: char
  INTEGER, PARAMETER :: lower_to_upper = IACHAR("A") - IACHAR("a")
  INTEGER :: i
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4 $//$Date: 2003/10/15 21:58:56 $'
  !     ------------------------------------------------------------------
  DO i=1,LEN(string)
     char = string(i:i)
     IF("A" <= char .AND. char <= "Z") THEN
        outstring(i:i) = char
     ELSEIF("a" <= char .AND. char <= "z") THEN
        outstring(i:i) = ACHAR(IACHAR(char) + lower_to_upper)
     ELSE      ! ... not alphabetic
        outstring(i:i) = char
     ENDIF
  END DO
END FUNCTION uppercase
