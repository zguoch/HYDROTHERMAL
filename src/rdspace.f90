SUBROUTINE rdspace(dumarray, array, fmt, utop, nl, kodd, convfac)
  ! ... Purpose:  To read input arrays for X, Y, and Z-spacing.
  ! ...           ZSPA (z-spacing) is either read from top to bottom or bottom
  ! ...           to top, depending upon UTOP
  USE machine_constants, ONLY: kdp
  USE f_units
  USE control, ONLY: ierr
  USE mesh, ONLY: nxyz
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nl
  REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: dumarray
  CHARACTER(LEN=4), INTENT(IN) :: array
  CHARACTER(LEN=80), INTENT(IN) :: fmt
  CHARACTER(LEN=3), INTENT(IN) :: utop
  INTEGER, INTENT(OUT) :: kodd
  REAL(KIND=kdp), INTENT(IN) :: convfac
  ! ...   DUMARRAY - array that is to be read (DX,DY, or DZ)
  ! ...   ARRAY - array name, either XSPA, YSPA, or ZSPA
  ! ...   FMT - specified format or method used to read array
  ! ...   UTOP  - indicates for FREE format if values are read from
  ! ...       top to bottom or bottom to top
  ! ...   NL  - Number of elements in array - either NX,NY, or NZ
  ! ...   KODD  =1 - array input as variable value array
  ! ...         =2 - array input as constant value array
  ! ...   CONVFAC - Conversion factor for array, from UNITFAC()
  !
!!$  CHARACTER(LEN=10), EXTERNAL :: uppercase
  CHARACTER(LEN=10) :: ustring
  INTEGER :: ll
  REAL(kind=kdp) :: dum1
  REAL(kind=kdp), DIMENSION(:), ALLOCATABLE :: fac
  INTEGER :: alloc_stat
  INTERFACE
     FUNCTION uppercase(string) RESULT(outstring)
       CHARACTER(LEN=*), INTENT(IN) :: string
       CHARACTER(LEN=LEN(string)) :: outstring
     END FUNCTION uppercase
  END INTERFACE
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.9 $//$Date: 2006/03/03 23:24:01 $'
  ALLOCATE(fac(nxyz), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: rdspace'
    ierr(181) = .TRUE.
    RETURN
  END IF
  !     ------------------------------------------------------------------
  ustring = ADJUSTL(fmt(1:10))
  ustring = uppercase(ustring(1:4))
  IF (ustring(1:4) == 'FREE') THEN          !  ... read array with FREE FORMAT
     kodd = 1
     IF (utop == 'TOP') THEN
        READ(fuins,*) (dumarray(ll), ll=nl,1,-1)
     ELSE
        READ(fuins,*) (dumarray(ll), ll=1,nl)
     END IF
  ELSE IF (ustring(1:4) == 'CONS') THEN     !  ... assign a CONSTANT VALUE to array
     kodd = 2
     READ(fuins,*) dum1
     DO  ll=1,nl
        dumarray(ll) = dum1
     END DO
  ELSE IF (ustring(1:4) == 'SCAL') THEN     !  ... free format with single scalar multiplier
     ! *** undocumented feature 
     kodd = 1
     READ(fuins,*) fac(1)
     IF (utop == 'TOP') THEN
        READ(fuins,*) (dumarray(ll), ll=nl,1,-1)
     ELSE
        READ(fuins,*) (dumarray(ll), ll=1,nl)
     END IF
     DO  ll = 1,nl
        dumarray(ll) = dumarray(ll)*fac(1)
     END DO
  ELSE IF (ustring(1:4) == 'FACT') THEN     !  ... free format with variable scalar multipliers
     ! *** undocumented feature 
     IF (utop == 'TOP') THEN
        READ(fuins,*) (dumarray(ll), fac(ll), ll=nl,1,-1)
     ELSE
        READ(fuins,*) (dumarray(ll), fac(ll), ll=1,nl)
     END IF
     DO  ll=1,nl
        dumarray(ll) = dumarray(ll)*fac(ll)
     END DO
  ELSE IF (ustring(1:4) == 'ROWC' .OR. ustring(1:4) == 'COLC' .OR.  &
       ustring(1:4) == 'NODE' .OR. ustring(1:4) == 'ROCK' .OR. ustring(1:4) == 'RTFT') THEN
     !  ... error if any of these formats are used
     !  ... constant value for each row or column, or individual node input
     !  ... format and error statements
     WRITE(fustdout,9005) '*** STOP - Error defining ',array,' ***',fmt(1:4),  &
          ' input is not valid for this parameter'
     WRITE(fuclog,9005) '*** STOP - Error defining ',array,' ***',fmt(1:4),  &
          ' input is not valid for this parameter'
9005 FORMAT (/a,a4,a/tr10,2a)
     ierr(25) = .TRUE.
  ELSE
     !  ... read the array with SPECIFIED FORMAT
     kodd = 1
     READ(fuins,fmt) (dumarray(ll),ll=1,nl)
  END IF
  !  ... check that cell dimensions are greater than zero
  DO  ll=1,nl
     IF (dumarray(ll) <= 0._kdp) THEN
        WRITE(fustdout,9010) array, ll, dumarray(ll)
        WRITE(fuclog,9010) array, ll, dumarray(ll)
9010    FORMAT (/, '*** STOP - Dimension error for ', a4, ' ***', /,  &
             10X, 'Node', i4, ' = ', 1P, g10.3)
        ierr(26) = .TRUE.
     END IF
  END DO
  ! ...  convert input units to internal cgs units (cm)
  DO  ll = 1, nl
     dumarray(ll) = dumarray(ll)*convfac
  END DO
  DEALLOCATE(fac)
END SUBROUTINE rdspace
