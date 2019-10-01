SUBROUTINE terminate_ht(flag)
  ! ... Terminates the simulation run invoking normal shut-down procedures
  ! ...      or error processing as necessary
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: flag
  !
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.6 $//$Date: 2003/10/15 22:08:57 $'
  !     ------------------------------------------------------------------
  !...
  ! ... Print current data arrays if an error exit
  IF(flag < 99) CALL pdata(flag)
  ! ... Normal closure steps
  CALL endsummary
  CALL closef(flag)
END SUBROUTINE terminate_ht
