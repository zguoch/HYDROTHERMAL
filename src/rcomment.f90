SUBROUTINE rcomment
  !     Purpose:  Reads the input file, strips out the comment lines, and
  !        writes the data to a temporary input file
  !        This temporary file is then read by the program
  USE f_units, ONLY: fuinc, fuins
  IMPLICIT NONE
  CHARACTER(10) :: buffer
  CHARACTER(len=:), ALLOCATABLE :: aline 
  INTEGER :: record_length
  INTEGER :: string_size
  !     ------------------------------------------------------------------------
  !
  DO
    ! Determine the length of a record (a line of theinput data file) by
    ! reading 10 characters at a time (without advancing to next record)
    ! until end of record (EOR) is reached.
    record_length = 0
    DO    
      READ (fuinc, '(a)', ADVANCE='NO', SIZE=string_size, EOR=10,   &
            END=20) buffer
      record_length = record_length + 10
    END DO
    10 record_length = record_length + string_size
    ! If the record is not empty, allocate a character string and reread
    ! the record into the character string. 
    IF (record_length > 0) THEN
      ALLOCATE (CHARACTER(LEN=record_length) :: aline)
      BACKSPACE (fuinc)
      READ (fuinc, '(a)') aline
      ! If the line is not a comment line, write it to the temporary
      ! input file
      IF (aline(1:1) /= '#') THEN
        WRITE (fuins, '(A)') aline
      END IF
      DEALLOCATE (aline)
    END IF
  END DO
  20 REWIND(fuins)
END SUBROUTINE rcomment
