SUBROUTINE wrnew(funit,dumarray,kodd,convfac)
  ! ... Purpose:  To write arrays.  Values in printout are changed from the
  ! ...           base units used in HYDROTHERM to match the user-specified
  ! ...           units
  ! ... Method:  After converting the units, the old wr subroutine is called.
  USE machine_constants, ONLY: kdp
  USE f_units
  USE control, ONLY: ierr
  USE mesh
  IMPLICIT NONE
  ! ...   FUNIT (input) - unit number to direct output
  ! ...   DUMARRAY (input) - array to be written
  ! ...   KODD (input) =1 - array input as variable value array
  ! ...         =2 - array input as uniform value array
  ! ...         =3 - array input as values uniform for each row
  ! ...         =4 - array input as values uniform for each column
  ! ...         =5 - array input as rocktype values
  ! ...         =50 - special code for writing Explorer plotfiles
  ! ...         =52 - special code for writing gnuplot plotfiles - unformatted
  ! ...         =53 - special code for writing gnuplot plotfiles - formatted
  ! ...         =54 - special code for writing idl plotfiles - unformatted
  ! ...         =55 - special code for writing idl plotfiles - formatted
  ! ...         =60 - special code for hydrostatic Pressure gradient
  ! ...         =70 - special code indicating P & H of all nodes are on the sat'd
  ! ...                  water or steam curve (hydrostatic P gradient for boiling
  ! ...                  point density)
  ! ...         =75 - special code indicating P & H of some columns are on the
  ! ...                  sat'd water or steam curve (hydrostatic P gradient for
  ! ...                  boiling point density)
  ! ...         =77 - special code indicating P & H of a node is on the sat'd
  ! ...                 water or steam curve (Note: no hydrostatic P gradient
  ! ...                 for boiling is implied)
  ! ...         =90 - special code indicating Y and/or Z permeabilities
  ! ...                 are identical to X permeabilities
  ! ...   CONVFAC (input) - Conversion factor from UNITFAC() for array
  INTEGER, INTENT(IN) :: funit
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: dumarray
  INTEGER, INTENT(IN) :: kodd
  REAL(KIND=kdp), INTENT(IN OUT) :: convfac
  !
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: dumarray2
  INTEGER :: alloc_stat
  INTERFACE
     SUBROUTINE wr(funit, dumarray, ikodd)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: funit
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: dumarray
       INTEGER, INTENT(IN) :: ikodd
     END SUBROUTINE wr
  END INTERFACE
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4 $//$Date: 2006/10/03 23:43:50 $'
  ALLOCATE(dumarray2(SIZE(dumarray)), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: wrnew'
    ierr(176) = .TRUE.
    RETURN
  END IF
  !     ------------------------------------------------------------------
  ! ... Convert from the base units (cgs) to the specified output units
  IF (convfac <= 0._kdp) THEN
     WRITE (fustdout,'(A)') 'Error in WRNEW: Conversion factor = 0'
     ierr(121) = .true.
     convfac = 1._kdp
  END IF
  dumarray2 = dumarray/convfac
  CALL wr(funit,dumarray2,kodd)
  DEALLOCATE(dumarray2)
END SUBROUTINE wrnew
