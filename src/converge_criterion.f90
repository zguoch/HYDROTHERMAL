SUBROUTINE converge_criterion(lconvrg, aconvrg, uconvrg)
  !     Purpose:  keep track of which criteria controlled NR convergence
!!$  USE f_units
  USE bc, ONLY: unconfined
  USE control
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: lconvrg
  CHARACTER (LEN=1), INTENT(OUT) :: aconvrg
  CHARACTER (LEN=13), INTENT(OUT) :: uconvrg
  !       LCONVRG - criteria number which met convergence specifications
  !              new convergence criteria, 1 of 3 possibilities
  !          1  if residuals are sufficienly small
  !          2  on/after 5 NRI, if absolute change in P&H is sufficiently small
  !          3  on/after 4 NRI, if % change in P&H is sufficiently small
  !               and residuals are within a factor of 100 of those on line 1
  !               and there has been a phase change
  !       ACONVRG  - single character indicating which criterion was met
  !       UCONVRG  - character string indicating which criterion was met
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.2.1.1 $//$Date: 2004/10/08 18:44:15 $'
  !     ------------------------------------------------------------------
  !...
  SELECT CASE(lconvrg)
  CASE (1)
     ! ... Convergence criteria met by mass and energy residuals (input parameter)
     iconvrg(1) = iconvrg(1) + 1
     aconvrg = 'R'
     uconvrg = 'Residuals'
  CASE (2)
     ! ... Convergence criteria met by absolute change in P&H
     ! ...      (hard coded in HYDROTHERM)
     iconvrg(2) = iconvrg(2) + 1
     aconvrg = 'C'
     uconvrg = 'AbsChgP&H'
  CASE (3)
     ! ... Convergence criteria met by % change in P&H (input parameter)
     iconvrg(3) = iconvrg(3) + 1
     IF(unconfined) THEN
        aconvrg = '%'
        uconvrg = 'AbsChgP,%ChgH'
     ELSE
        aconvrg = '%'
        uconvrg = '%Chg P&H'
     END IF
  END SELECT
END SUBROUTINE converge_criterion
