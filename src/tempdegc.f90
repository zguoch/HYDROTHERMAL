SUBROUTINE tempdegc(p,h,tc)
  ! ... Purpose:  To compute temperature from pressure and enthalpy
  ! ... Note:  Index of thermodynamic phase region for nodes must be current before
  ! ...      running this subroutine
  USE machine_constants, ONLY: kdp
  USE f_units
  USE mesh
  USE variables, ONLY: ind
  IMPLICIT NONE
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
  REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: tc
  !
  REAL(KIND=kdp) :: dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8,  &
       dum9, duma, dumb, dumc, dumd, dume
  INTEGER :: ic, icdum
  LOGICAL :: errflag
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.4 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  DO  icdum = 1,npiccons
     ic = npic(icdum)
     IF (ind(ic) == 2) THEN           ! ... for two phase
        CALL phsbndy3(p(ic), dum1, dum2, dum3, dum4, dum5, dum6, dum7,  &
             dum8, dum9, duma, dumb, dumc, tc(ic), dumd, dume)
     ELSE          ! ... for water, super-heated steam, supercritcal fluid, air-water
        CALL tblookup(p(ic), h(ic), dum1, dum2, dum3, tc(ic), dum4,  &
             dum5, dum6, dum7, dum8, errflag)
        IF(errflag) THEN
           WRITE (fustdout,9005)  &
                '***** ERROR: Pressure-Enthalpy out of table range for cell no. ',ic,' *****'
           WRITE (fuclog,9005)  &
                '***** ERROR: Pressure-Enthalpy out of table range error for cell no. ',ic,' *****'
9005       FORMAT(a,i6,a)
        END IF
     END IF
  END DO
END SUBROUTINE tempdegc

