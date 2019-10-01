SUBROUTINE plotexpl
  ! ... Purpose:  To write results for the current time-level on disk for
  ! ...       later visualization using the SGI Explorer program
  !
  ! ... Method:  For each time step that is to be plotted, several files
  ! ...       are written, each with a number designation at the end
  ! ...       of the file name to indicate the relative time step.
  ! ...       Separate files are written for all the variable that are to
  ! ...       be plotted.  Each time step, this subroutine produces
  ! ...       separate files for pressure, enthalpy, temperature,
  ! ...       water saturation.  Other files are written for plotting
  ! ...       water and steam velocity vectors; these files include
  ! ...       X, Y, and Z velocities for water and steam.
  ! ... note: IOPTPR(1) = 1 for Explorer plotfiles
  ! ... note: nodes that are in the supercritical region (P>Pc) are
  ! ...       assigned a sat'n of -1 in this subroutine for ploting purposes.
  USE machine_constants, ONLY: kdp
  USE f_units
  USE control
  USE mesh
  USE variables
  USE units
  IMPLICIT NONE
  CHARACTER (LEN=80) :: ufile
  INTEGER :: ic, icxf, icxl, j, k, lchar1, lchar2, lchar3
  ! ...   ipltkodd = 52 for binary plotfile  = 53 for ASCII plotfile;
  LOGICAL :: lprint
  REAL(kind=kdp), DIMENSION(:), ALLOCATABLE :: satnw2
  INTEGER :: alloc_stat
  INTERFACE
     SUBROUTINE wr(funit, dumarray, ikodd)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: funit
       REAL(kind=kdp), DIMENSION(:), INTENT(IN) :: dumarray
       INTEGER, INTENT(IN) :: ikodd
     END SUBROUTINE wr
  END INTERFACE
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.2 $//$Date: 2007/12/20 20:32:47 $'
  ALLOCATE(satnw2(nxyz), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: plotexpl'
    ierr(186) = .TRUE.
    RETURN
  END IF
  !     ------------------------------------------------------------------
  !...
!!$  !---set counter for "relative time step" for filename suffix
!!$  lpltexpl = lpltexpl + 1
!!$  IF (lpltexpl >= 100) THEN
!!$     IF (lpltexpl <= 110) WRITE (fustdout, '(5X,A)')  &
!!$          '**** Limit of 100 plotfiles exceeded ****'
!!$     WRITE(fuclog, '(10X,A)')  &
!!$          'Limit of 100 plotfiles exceeded, increase limit in plotexpl.f'
!!$     RETURN
!!$  END IF
!!$  !---convert value of LPLTEXPL to a character
!!$  lchar1 = IDINT(DBLE(lpltexpl)/100.0_KDP)
!!$  lchar2 = IDINT(DBLE(lpltexpl-lchar1*100)/10.0_KDP)
!!$  lchar3 = lpltexpl - (lchar1*100+lchar2*10)
!!$  ufile = CHAR(lchar1+48)//CHAR(lchar2+48)//CHAR(lchar3+48)
!!$  ! ...    open files
!!$  CALL gfiles(6, ufile)
  ! ... files opened in gfiles
  IF(ABS(print_plotscalar) > 0._kdp) THEN
     CALL print_control(print_plotscalar,time,ltime,tchg,timprpscal,lprint)
     IF(lprint) THEN
        !---print arrays
        !----- temperature -----
        REWIND(50)
        CALL plotheader(50)
        CALL wr(50,tc,50)
!!$  CLOSE(50)
        !----- pressure -----
!!$  !---create a pressure array in atm instead of dynes for plotting
!!$  !     DO 5 ic = 1, Nxyz
!!$  !       p2(ic) = atm(P(ic))
!!$  !   5 CONTINUE
!!$  REWIND(51)
        CALL plotheader(51)
        CALL wr(51,p,50)
!!$  CLOSE (51)
        !----- water saturation -----
        !---create a modified sat'n array SATNW2 so that all supercrital
        !---nodes (P>Pc) have a saturation of -1 instead of 0.5
        DO  ic = 1, nxyz
           IF (p(ic) > 220.55e6_kdp) THEN
              !              use the following if supercritical includes T>Tc
              !       IF ( (P(IC).GT.220.55D6.AND.SATNW(IC).EQ.0.5_KDP).OR.
              !    &       (TC(IC).GE.374.0_KDP.AND.SATNW(IC).EQ.0.0_KDP) ) THEN
              satnw2(ic) = -1.0_KDP
           ELSE
              satnw2(ic) = satnw(ic)
           END IF
        END DO
!!$  REWIND(52)
        CALL plotheader(52)
        CALL wr(52,satnw2,50)
!!$  CLOSE(52)
        !----- enthalpy -----
!!$  REWIND (53)
        CALL plotheader(53)
        CALL wr(53,h,50)
!!$  CLOSE(53)
     END IF
  END IF
  IF(ABS(print_plotvector) > 0._kdp) THEN
     CALL print_control(print_plotvector,time,ltime,tchg,timprpvec,lprint)
     IF(lprint) THEN
        !----- velocity of water and steam vectors -----
        !---write file (unit=54)with (X,Z) pairs for water velocity vectors
        !---write file (unit=55)with (X,Z) pairs for steam velocity vectors
!!$  REWIND(54)
        CALL plotheader(54)
!!$  REWIND (55)
        CALL plotheader(55)
        DO  j=1,ny
           DO  k=1,nz
              icxf = (k-1)*nx + 1 + (j-1)*nx*nz
              icxl = icxf + nx - 1
              WRITE(54, 9005) (xwvel(ic), zwvel(ic), ic=icxf, icxl)
              WRITE(55, 9005) (xsvel(ic), zsvel(ic), ic=icxf, icxl)
9005          FORMAT (1P, 100E11.3)
           END DO
        END DO
!!$  CLOSE(54)
!!$  CLOSE(55)
     END IF
  END IF
  DEALLOCATE(satnw2)

CONTAINS
  SUBROUTINE plotheader(lunit)
    !---header information for plot files
    USE mesh
    USE units
    USE variables
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: lunit
    !
    INTEGER :: i, j, k
    !     ------------------------------------------------------------------
    !...
    WRITE (lunit, '(3I5)') nx, nz, ny
    WRITE (lunit, 9005) (xx(i), i=1, nx)
    WRITE (lunit, 9005) (zz(k), k=1, nz)
    WRITE (lunit, 9005) (yy(j), j=1, ny)
    WRITE (lunit, 9005)
9005 FORMAT (1P, 50E11.3, (50E11.3))
    WRITE (LUNIT,'(1P,G12.6)') year(TIME)
  END SUBROUTINE plotheader

END SUBROUTINE plotexpl

