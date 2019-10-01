SUBROUTINE wr(funit, dumarray, ikodd)
  !     Purpose:  To write arrays.  Values in printout are changed from the
  !               base units used in HYDROTHERM to match the user-specified
  !               units
  USE machine_constants, ONLY: kdp
  USE control
  USE mesh
  IMPLICIT NONE
  !       FUNIT  - unit number to direct output
  !       DUMARRAY  - array to be written
  !       IKODD  =1 - array input as variable value array
  !             =2 - array input as uniform value array
  !             =3 - array input as values uniform for each row
  !             =4 - array input as values uniform for each column
  !             =5 - array input as rocktype values
  !             =50 - special code for writing Explorer plotfiles
  !             =52 - special code for writing gnuplot plotfiles - unformatted
  !             =53 - special code for writing gnuplot plotfiles - formatted
  !             =60 - special code for hydrostatic Pressure gradient
  !             =70 - special code indicating P & H of all nodes are on the sat'd
  !                      water or steam curve (hydrostatic P gradient for boiling
  !                      point density)
  !             =75 - special code indicating P & H of some columns are on the
  !                      sat'd water or steam curve (hydrostatic P gradient for
  !                      boiling point density)
  !             =77 - special code indicating P & H of a node is on the sat'd
  !                     water or steam curve (Note: no hydrostatic P gradient
  !                     for boiling is implied)
  !             =90 - special code indicating Y and/or Z permeabilities
  !                     are identical to X permeabilities
  !       CONVFAC  - Conversion factor from UNITFAC() for array
  INTEGER, INTENT(IN) :: funit
  REAL(kind=kdp), DIMENSION(:), INTENT(IN) :: dumarray
  INTEGER, INTENT(IN) :: ikodd
  !
  CHARACTER(LEN=80) :: cfmt
  CHARACTER(LEN=42), PARAMETER :: cfmt1="('Row',i3,' -> ',1P,6G12.4)             ",  &
      cfmt2="('++ Row',i3,' ++'/1P,7G12.4)           ",  &
      cfmt3="('++ Row',i3,' ++'/1P,7G12.4/(7G12.4))  ",  &
      cfmt4="('Row',i3,' -> ',1P,11G12.4)            ",  &
      cfmt5="('++ Row',i3,' ++'/1P,12G12.4)          ",  &
      cfmt6="('++ Row',i3,' ++'/1P,12G12.4/(12G12.4))",  &
      cfmt7="('Row',i3,' -> ',1P,13G12.4)            ",  &
      cfmt8="('++ Row',i3,' ++'/1P,14G12.4)          ",  &
      cfmt9="('++ Row',i3,' ++'/1P,14G12.4/(14G12.4))"
  INTEGER :: i, ic, icxf, icxl, iczf, j, k, kodd
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.7 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  kodd = ikodd
  IF (kodd == 0) kodd = 1
  IF (kodd == 50) THEN
     !                   ***** write array for explorer plotfile *****
     DO  j = 1, ny
        DO  k = 1, nz
           icxf = (k-1)*nx + 1 + (j-1)*nxz
           icxl = icxf + nx - 1
           WRITE(funit,9035) (dumarray(ic),ic=icxf,icxl)
        END DO
     END DO
  ELSE IF (kodd == 52) THEN     ! ... array for gnuplot plotfile - unformatted
     DO  j = 1,ny
        DO  k = 1,nz
           icxf = (k-1)*nx + 1 + (j-1)*nxz
           icxl = icxf + nx - 1
           WRITE(funit) (dumarray(ic),ic=icxf,icxl)
        END DO
     END DO
  ELSE IF (kodd == 53) THEN     ! ... array for gnuplot plotfile - formatted
     DO  j = 1,ny
        DO  k = 1,nz
           icxf = (k-1)*nx + 1 + (j-1)*nxz
           icxl = icxf + nx - 1
           WRITE(funit,9040) (dumarray(ic),ic=icxf,icxl)
        END DO
     END DO
  ELSE IF (kodd == 54) THEN     ! ... array for idl plotfile - unformatted
     DO  k = 1,nz
        DO  j = 1,ny
           icxf = (k-1)*nx + 1 + (j-1)*nxz
           icxl = icxf + nx - 1
           WRITE(funit) (dumarray(ic),ic=icxf,icxl)
        END DO
     END DO
  ELSE IF (kodd == 55) THEN     ! ... array for idl plotfile - formatted
     DO  k = 1,nz
        DO  j = 1,ny
           icxf = (k-1)*nx + 1 + (j-1)*nxz
           icxl = icxf + nx - 1
           WRITE(funit,9040) (dumarray(ic),ic=icxf,icxl)
        END DO
     END DO
  ELSE IF (kodd < 2 .OR. kodd > 4) THEN
     !.....array with variable values; KODD= 5 or 0 or 1
     !             format for narrow (80 col) print option
     IF (ioptpr(2) == 0) THEN
        IF (nx < 7) THEN
           cfmt = cfmt1
        ELSE IF (nx == 7) THEN
           cfmt = cfmt2
        ELSE
           cfmt = cfmt3
        END IF
        !             format for wide (132 col) print option
     ELSE IF (ioptpr(2) == 1) THEN
        IF (nx < 12) THEN
           cfmt = cfmt4
        ELSE IF (nx == 12) THEN
           cfmt = cfmt5
        ELSE
           cfmt = cfmt6
        END IF
        !             format for ex-wide (159 col) print option
     ELSE IF (ioptpr(2) == 2) THEN
        IF (nx < 14) THEN
           cfmt = cfmt7
        ELSE IF (nx == 14) THEN
           cfmt = cfmt8
        ELSE
           cfmt = cfmt9
        END IF
     END IF
     DO  j = 1, ny
        WRITE (funit, '(1X)')
        IF (ny > 1) WRITE (funit, 9010) j
        DO  k = nz, 1, -1
           icxf = (k-1)*nx + 1 + (j-1)*nxz
           icxl = icxf + nx - 1
           WRITE(funit,cfmt) k,(dumarray(ic),ic=icxf,icxl)
        END DO
     END DO
     !.....New array printing routine to be installed later
  ELSE IF (kodd == 2) THEN
     !.....uniform value array;KODD= 2
     WRITE(funit,9005) dumarray(1)
9005 FORMAT (10X, 'Uniform Value:  ', 1P,g11.4)
  ELSE IF (kodd == 3) THEN
     !.....array with rows of uniform value; KODD= 3
     DO  j = 1, ny
        WRITE (funit, 9015)
9015    FORMAT (15X, 'Each Row Has Uniform Values')
        IF (ny > 1) WRITE (funit, 9010) j
        DO  k = nz, 1, -1
           icxf = (k-1)*nx + 1 + (j-1)*nxz
           WRITE(funit,9020) k,dumarray(icxf)
        END DO
     END DO
  ELSE IF (kodd == 4) THEN
     !.....array with columns of uniform value; KODD= 4
     DO  j = 1, ny
        WRITE (funit, 9025)
9025    FORMAT (15X, 'Each Column has Uniform Values')
        IF (ny > 1) WRITE (funit, 9010) j
        DO  i = 1, nx
           iczf = i + (j-1)*nxz
           WRITE(funit,9030) i,dumarray(iczf)
        END DO
     END DO
  END IF

  9010 FORMAT (10X, '--- Slice', i4, ' ---')
  9020 FORMAT ('Row', i4, ' -> ', 1P, g11.4)
  9030 FORMAT ('Column', i4, ' -> ', 1P, g11.4)
  9035 FORMAT (10(1pe11.3))
  9040 FORMAT (14(1pe11.3))
END SUBROUTINE wr
