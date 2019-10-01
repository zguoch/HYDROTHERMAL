SUBROUTINE rd(dumarray, fmt, utop, kodparm, kodd, dumarray2)
  ! ... Purpose:  To read input arrays one slice at a time, row by row,
  ! ...     from bottom to top or top to bottom
  USE machine_constants, ONLY: kdp
  USE f_units
  USE bc, ONLY: ibc, prflag
  USE control
  USE mesh
  USE parameters
  USE units
  IMPLICIT NONE
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: dumarray
  CHARACTER(LEN=80), INTENT(IN) :: fmt
  CHARACTER(LEN=3), INTENT(IN) :: utop
  INTEGER, INTENT(IN) :: kodparm
  INTEGER, INTENT(OUT) :: kodd
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT), OPTIONAL :: dumarray2
  ! ...   DUMARRAY  - array that is to be read
  ! ...   FMT  - specified format or method used to read array
  ! ...   UTOP - indicates if values are read from top to bottom or bottom to top
  ! ...   KODPARM - Parameter index
  ! ...         1 - Porosity
  ! ...         2 - X permeability
  ! ...         3 - Y permeability
  ! ...         4 - Z permeability
  ! ...         5 - Thermal Conductivity
  ! ...         6 - Specific Heat
  ! ...         7 - Rock Density
  ! ...         8 - Rock Compressibility
  ! ...        NOTE: 9 - Conductive Heat input is not input with this subroutine
  ! ...        10 - Enthalpy
  ! ...        11 - Pressure
  ! ...        15 - Temperature
  ! ...   KODD  =1 - array input as variable value array
  ! ...         =2 - array input as constant value array
  ! ...         =3 - array input as values constant for each row
  ! ...         =4 - array input as values constant for each column
  ! ...          5 - array input as rocktype values
  ! ...          6 - array input as rocktype with function of temperature
  ! ...          7 - array input as rocktype with function of time
  ! ...         =50 - code for writing Explorer plotfiles
  ! ...         =60 - code for hydrostatic Pressure gradient
  ! ...         =70 - code indicating P & H of all nodes are on the sat'd
  ! ...                  water or steam curve (hydrostatic P gradient for boiling
  ! ...                  point density)
  ! ...         =75 - code indicating P & H of some columns are on the
  ! ...                  sat'd water or steam curve (hydrostatic P gradient for
  ! ...                  boiling point density)
  ! ...         =77 - code indicating P & H of a node is on the sat'd
  ! ...                 water or steam curve (Note: no hydrostatic P gradient
  ! ...                 for boiling is implied)
  ! ...         =90 - code indicating Y and/or Z permeabilities
  ! ...                 vary by some factor relative to X permeabilities
  !
  CHARACTER(LEN=4), DIMENSION(mprxparm), PARAMETER :: chrarry=(/'PORO', 'XPER', 'YPER', 'ZPER', &
       'THER', 'SPEC', 'DENS', 'COMP'/)
!!$  CHARACTER(LEN=10), EXTERNAL :: uppercase
  CHARACTER(LEN=10) :: ustring
  INTEGER :: i, ic, icdum, icxf, icxl, idum1, idum2, ii, iii, ik, irx, j,  &
       k, kk, lll
  REAL(KIND=kdp) :: dum1
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: dumnx
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: dumnz
  INTEGER :: alloc_stat
  INTERFACE
     FUNCTION uppercase(string) RESULT(outstring)
       CHARACTER(LEN=*), INTENT(IN) :: string
       CHARACTER(LEN=LEN(string)) :: outstring
     END FUNCTION uppercase
  END INTERFACE
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.19 $//$Date: 2007/08/03 23:42:47 $'
  ALLOCATE(dumnx(nx), dumnz(nz), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: rd'
    ierr(182) = .TRUE.
    RETURN
  END IF
  !     ------------------------------------------------------------------
  ! ...  use ROCK TYPE input for arrays 
  ustring = ADJUSTL(fmt(1:10))
  ustring = uppercase(ustring(1:10))
  IF (ustring(1:4) == 'ROCK') THEN          ! ... Input by rock type or zone
     IF (kodparm < 1 .OR. kodparm > 8) THEN
        ! ... Error if ROCK type input is not available for current parameter
        WRITE (fustdout,9015) fmt(1:4),kodparm
        WRITE (fuclog,9015) fmt(1:4),kodparm
9015    FORMAT (/ '*** ERROR - Problem with array input ***' /10X,a,  &
             ' input is not valid for parameter, kodparm= ',i2)
        ierr(21) = .TRUE.
     END IF
     kodd = 5
     DO  lll = 1, nrxtype
        READ(fuins,*,ERR=200) irx, dum1
        IF (ilrxtype(irx)) THEN      ! ... convert from user units to internal units (cgs)
                                     ! ...      and assign value to rock type
           rxparm(kodparm,irx) = dum1*unitfac(kodparm)
        ELSE                ! ... invalid rock type - write error statement
           WRITE (fustdout,9010) chrarry(kodparm)
           WRITE (fustdout,9020) fmt(1:7)
9020       FORMAT (tr10, 'Problem with ', a, ' input')
           WRITE (fustdout,9030) irx
           WRITE (fuclog,9010) chrarry(kodparm)
           WRITE (fuclog,9020) fmt(1:7)
           WRITE (fuclog,9030) irx
9030       FORMAT (tr10,'Rock type ',i2,' not used in this simulation')
           ierr(24) = .TRUE.
        END IF
     END DO
     ! ... Load values for the rock types into array
     DO  icdum=1,npiccons
        ic = npic(icdum)
        irx = icrxtype(ic)
        dumarray(ic) = rxparm(kodparm,irx)
     END DO
  ELSEIF ((informt >= 3 .AND. ustring(1:7) == 'RTFTEMP') .OR.    &
          (informt <= 2 .AND. ustring(1:4) == 'RXFT')) THEN
     ! ... Variables for Rock Parameters that are Functions of Temperature
     IF (kodparm < 1 .OR. kodparm > 8) THEN      ! ... parameter not appropriate for RTFTEMP input
        WRITE (fustdout,9015) fmt(1:7),kodparm
        WRITE (fuclog,9015) fmt(1:7),kodparm
        ierr(22) = .TRUE.      ! ... Error; ROCK type input not available for this parameter
     END IF
     kodd = 6
     DO  lll = 1,nrxtype
        READ(fuins,*,ERR=200) irx, idum1, idum2
        ! ... IDUM1 is the option number
        ! ... IDUM2 is the number of points defining f(Temperature)
        IF (ilrxtype(irx)) THEN
           irxftopt(kodparm,irx,1) = idum1
           irxftopt(kodparm,irx,2) = idum2
           ! ... check that option number is valid
           IF (idum1 < 1 .OR. idum1 > 4) THEN
              WRITE (fustdout,9010) chrarry(kodparm)
              WRITE (fustdout,'(tr10,A,I2)') 'Invalid option number: ',idum1
              WRITE (fuclog,9010) chrarry(kodparm)
              WRITE (fuclog,'(tr10,A,I2)') 'Invalid option number: ',idum1
              ierr(61) = .TRUE.
           END IF
           ! ...                 check that number of points is > 1 or <= 4 
           IF (idum2 > 4) THEN
              WRITE (fustdout,9010) chrarry(kodparm)
              WRITE (fustdout,'(A)')  &
                   'Number of points defining RTFTEMP function greater than maximum allowed (>4)'
              WRITE (fuclog,9010) chrarry(kodparm)
              WRITE (fuclog,'(A)')  &
                   'Number of points defining RTFTEMP function greater than maximum allowed (>4)'
              ierr(62) = .TRUE.
           END IF
           IF (idum1 >= 2 .AND. idum2 <= 1) THEN
              WRITE (fustdout,9010) chrarry(kodparm)
              WRITE (fustdout,'(A)')  &
                   'Number of points defining RTFTEMP function less than minimum required (<2)'
              WRITE (fuclog,9010) chrarry(kodparm)
              WRITE (fuclog,'(A)')  &
                   'Number of points defining RTFTEMP function less than minimum required (<2)'
              ierr(62) = .TRUE.
           END IF        
           IF (idum1 == 1) THEN    ! ... function #1 - parameter is uniform (not a function of Temperature)
              READ(fuins,*) dum1
              ! ... convert from user specified units to internal units (cgs)
              ! ...     and assign value to rock type
              rxparm(kodparm,irx) = dum1*unitfac(kodparm)
              rxftparm(kodparm,irx,1) = dum1*unitfac(kodparm)
              ! ... assign values to the array for function 1
              DO  icdum = 1, npiccons
                 ic = npic(icdum)
                 irx = icrxtype(ic)
                 dumarray(ic) = rxparm(kodparm,irx)      ! ... no need to calculate as f(T)
              END DO
           ELSE IF (idum1 == 2) THEN     ! ... function #2 - parameter is a step function of Temperature
              READ(fuins,*) (rxftparm(kodparm,irx,ii), ii=1, 2*idum2-1)
              ! ...            convert from user specified units to internal units (cgs)
              ! ...              NOTE: no conversion is done here for the Temperature in (deg C)
              DO  ii = 1,idum2
                 rxftparm(kodparm,irx,2*ii-1) = rxftparm(kodparm,irx,2*ii-1)*unitfac(kodparm)
              END DO
              rxparm(kodparm,irx) = -1.0_kdp
              ! ...               when 3 or more step functions are input check that
              ! ...                 values are input from low to high temperature.
              ! ...                 it is impossible to check with 2 steps because this
              ! ...                   only requires 1 Temperature value to define the 2 steps
              IF (idum2 >= 3) THEN
                 DO  ii = 1,idum2 - 2
                    IF (rxftparm(kodparm,irx,2*ii) >=  &
                         rxftparm(kodparm,irx,2*ii+2)) THEN
                       WRITE (fustdout,'(//,A,A)') 'ERROR: Check that parameter ',  &
                            'values are input in order of increasing '// 'temperature'
                       WRITE (fuclog,'(//,A,A)') 'ERROR: Check that parameter ',  &
                            'values are input in order of increasing '// 'temperature'
                       ierr(23) = .TRUE.
                    END IF
                 END DO
              END IF
           ELSE IF (idum1 == 3 .OR. idum1 == 4) THEN
              ! ... function #3 - log parameter is a linear function of Temperature
              ! ...     and function #4 - parameter is a linear function of Temperature
           ! *** disable function #3 for specific heat because computation in subroutine enthrock has 
           !     not been implemented ***
              IF (idum1 == 3 .AND. kodparm == 6) THEN
                  WRITE (fustdout,'(//,A,A)') 'ERROR: Function no. 3 for temperature-depencence of ',  &
                 'specific heat is currently disabled. Please use an alternative function.'
                  WRITE (fuclog,'(//,A,A)') 'ERROR: Function no. 3 for temperature-depencence of ',  &
                 'specific heat is currently disabled. Please use an alternative function.'
                  ierr(172) = .TRUE.
              END IF
              READ(fuins,*) (rxftparm(kodparm,irx,ii), ii=1, 2*idum2)
              ! ...            convert from user specified units to internal units (cgs)
              ! ...              NOTE: no conversion is done here for the Temperature (deg C)
              DO  kk = 1,idum2
                 rxftparm(kodparm,irx,2*kk-1) = rxftparm(kodparm,irx,2*kk-1)*unitfac(kodparm)
              END DO
              rxparm(kodparm,irx) = -1.0_kdp
              ! ... when 2 or more linear segments are input check that
              ! ...      values are input from low to high temperature
              IF (idum2 >= 2) THEN
                 DO  ii = 1,idum2 - 1
                    IF (rxftparm(kodparm,irx,2*ii) >=  &
                         rxftparm(kodparm,irx,2*ii+2)) THEN
                       WRITE (fustdout,'(//,A,A)') 'ERROR: Check that parameter ',  &
                            'values are input in order of increasing temperature'
                       WRITE (fuclog,'(//,A,A)') 'ERROR: Check that parameter ',  &
                            'values are input in order of increasing temperature'
                       ierr(23) = .TRUE.
                    END IF
                 END DO
              END IF
           END IF
           ! ...           stop if rock type is not defined - write error statements
        ELSE
           WRITE (fustdout,9010) chrarry(kodparm)
           WRITE (fustdout,9020) fmt(1:7)
           WRITE (fustdout,9030) irx
           WRITE (fuclog,9010) chrarry(kodparm)
           WRITE (fuclog,9020) fmt(1:7)
           WRITE (fuclog,9030) irx
           ierr(24) = .TRUE.
        END IF
     END DO
  ELSEIF (ustring(1:7) == 'RTFTIME') THEN
     ! ------- VARIABLES for ROCK PARAMETERS that are FUNCTIONS of TIME ----
     ! ...       stop if parameter is not appropriate for RTFTIME input
     IF (kodparm < 1 .OR. kodparm > 4) THEN      ! ... only permeability and porosity
        WRITE (fustdout,9015) fmt(1:7),kodparm
        WRITE (fuclog,9015) fmt(1:7),kodparm
        ierr(76) = .TRUE.
     END IF
     kodd = 7
     DO  lll = 1, nrxtype
        READ(fuins,*,ERR=200) irx, idum1, idum2
        ! ... IDUM1 is the option number
        ! ... IDUM2 is the number of points defining f(time)
        IF (ilrxtype(irx)) THEN
           irxftimopt(kodparm,irx,1) = idum1
           irxftimopt(kodparm,irx,2) = idum2
           ! ... check that option number is valid
           IF (idum1 < 1 .OR. idum1 > 2) THEN
              WRITE (fustdout,9010) chrarry(kodparm)
              WRITE (fustdout,'(10X,A,I2,A)') 'Function number ',idum1,  &
                   ' is not valid'
              WRITE (fuclog,9010) chrarry(kodparm)
              WRITE (fuclog,'(10X,A,I2,A)') 'Function number ',idum1,  &
                   ' is not valid'
              ierr(77) = .TRUE.
           END IF
           IF (idum2 > 4) THEN
              WRITE (fustdout,9010) chrarry(kodparm)
              WRITE (fustdout,'(A)')  &
                   'Number of points defining RTFTIME function greater than maximum allowed (>4)'
              WRITE (fuclog,9010) chrarry(kodparm)
              WRITE (fuclog,'(A)')  &
                   'Number of points defining RTFTIME function greater than maximum allowed (>4)'
              ierr(78) = .TRUE.
           END IF
           IF (idum1 >= 2 .AND. idum2 <= 1) THEN
              WRITE (fustdout,9010) chrarry(kodparm)
              WRITE (fustdout,'(A)')  &
                   'Number of points defining RTFTIME function less than minimum required (<2)'
              WRITE (fuclog,9010) chrarry(kodparm)
              WRITE (fuclog,'(A)')  &
                   'Number of points defining RTFTIME function less than minimum required (<2)'
              ierr(78) = .TRUE.
           END IF
           IF (idum1 == 1) THEN
              ! ... option 1 - constant parameter (not a function of time)
              READ(fuins,*) dum1
              ! ... convert from user specified units to internal units (cgs)
              ! ...     and assign value to rock type
              rxparm(kodparm,irx) = dum1*unitfac(kodparm)
              rxftimparm(kodparm,irx,1) = dum1*unitfac(kodparm)
              ! ... assign values to the array
              DO  icdum=1,npiccons
                 ic = npic(icdum)
                 irx = icrxtype(ic)
                 dumarray(ic) = rxparm(kodparm,irx)
              END DO
           ELSEIF (idum1 == 2) THEN
              ! ... option 2 - parameter is a linear function of time
              READ(fuins,*) (rxftimparm(kodparm,irx,ii), ii=1,2*idum2)
              ! ... convert from user specified units to internal units (cgs)
              ! ...      for parameter and time
              DO  ii=1,idum2
                 rxftimparm(kodparm,irx,2*ii-1) = rxftimparm(kodparm,irx,2*ii-1)*unitfac(kodparm)
                 IF(iyr == 2) rxftimparm(kodparm,irx,2*ii) = sec(rxftimparm(kodparm,irx,2*ii))
              END DO
              rxparm(kodparm,irx) = -1.0_kdp
              ! ... when 2 or more linear segments are input check that
              ! ...      values are input from low to high time
              IF (idum2 >= 2) THEN
                 DO  ii=1,idum2-1
                    IF(rxftimparm(kodparm,irx,2*ii) >= rxftimparm(kodparm,irx,2*ii+2)) THEN
                       WRITE(fustdout,'(//A,A)') 'ERROR: Check that parameter ',  &
                            'values are input in order of increasing time'
                       WRITE(fuclog,'(//A,A)') 'ERROR: Check that parameter ',  &
                            'values are input in order of increasing time'
                       ierr(79) = .TRUE.
                    END IF
                 END DO
              END IF
           END IF
           ! ...           stop if rock type is not defined - write error statements
        ELSE
           WRITE (fustdout,9010) chrarry(kodparm)
           WRITE (fustdout,9020) fmt(1:7)
           WRITE (fustdout,9030) irx
           WRITE (fuclog,9010) chrarry(kodparm)
           WRITE (fuclog,9020) fmt(1:7)
           WRITE (fuclog,9030) irx
           ierr(24) = .TRUE.
        END IF
     END DO
  ELSEIF (ustring(1:4) == 'FREE') THEN
     ! ...   ----------------- read array values with FREE FORMAT ------------------
     kodd = 1
     IF (utop == 'TOP') THEN
        DO  j = 1,ny
           DO  k = nz,1,-1
              icxf = (k-1)*nx + 1 + (j-1)*nxz
              icxl = icxf + nx - 1
              READ(fuins,*) (dumarray(ic), ic=icxf,icxl)
           END DO
        END DO
     ELSE
        READ(fuins,*) (dumarray(ic), ic=1,nxyz)
     END IF
     ! ...            convert from user specified units to internal units (cgs)
     DO ic = 1,nxyz
        dumarray(ic) = dumarray(ic)*unitfac(kodparm)
     END DO
  ELSEIF (ustring(1:4) == 'CONS') THEN
     ! ...   ----------------- assign a CONSTANT VALUE to array  ------------------
     kodd = 2
     READ(fuins, *) dum1
     ! ...            convert from user specified units to internal units (cgs)
     ! ...              and assign value to entire array
     dum1 = dum1*unitfac(kodparm)
     DO  ic = 1,nxyz
        dumarray(ic) = dum1
     END DO
     ! ...                    set values for rock types, if applicable
     IF(kodparm <= 8) THEN
        DO  lll = 1,nrxtype
           irx = irxused(lll)
           rxparm(kodparm,irx) = dum1
        END DO
     END IF
  ELSEIF (ustring(1:4) == 'ROWC') THEN
     ! ...   ------------ assign a CONSTANT VALUE to each ROW in array  -----------
     kodd = 3
     IF (utop == 'TOP') THEN
        READ(fuins,*) (dumnz(k), k=nz,1,-1)
     ELSE
        READ(fuins,*) (dumnz(k), k=1,nz)
     END IF
     DO  j=1,ny
        DO  k=1,nz
           icxf = (k-1)*nx + 1 + (j-1)*nxz
           icxl = icxf + nx - 1
           DO  ic=icxf,icxl
              dumarray(ic) = dumnz(k)
           END DO
        END DO
     END DO
     ! ...            convert from user specified units to internal units (cgs)
     DO ic = 1,nxyz
        dumarray(ic) = dumarray(ic)*unitfac(kodparm)
     END DO
  ELSEIF (ustring(1:4) == 'COLC') THEN
     ! ...   ------------ assign a CONSTANT VALUE to each COLUMN in array  --------
     kodd = 4
     READ(fuins,*) (dumnx(i), i=1,nx)
     DO  j = 1,ny
        DO  i = 1,nx
           DO  k = 1,nz
              ic = (k-1)*nx + i + (j-1)*nxz
              dumarray(ic) = dumnx(i)
           END DO
        END DO
     END DO
     ! ...            convert from user specified units to internal units (cgs)
     DO ic = 1,nxyz
        dumarray(ic) = dumarray(ic)*unitfac(kodparm)
     END DO
  ELSEIF (ustring(1:4) == 'NODE') THEN
     ! ...   ------------ read values for a subset of selected nodes ---------
     kodd = 1
     prflag = -1                     ! ... turn off node value printing
     DO 
        READ(fuins,*) i, j, k, dum1
        IF (i == 0 .AND. j == 0 .AND. k == 0) EXIT     ! ... terminate input
        IF (dum1 == -77._kdp) kodd = 77
        ic = (k-1)*nx + i + (j-1)*nxz
        ik = (k-1)*nx + i
        ! ...           convert from user specified units to internal units (cgs)
        ! ...             and assign value to node
        dumarray(ic) = dum1*unitfac(kodparm)
        IF(PRESENT(dumarray2)) dumarray2(ic) = -99._kdp     ! ... mark the unused input variable
        prflag(ik,j) = 1              ! ... turn on printing for this node
     END DO
  ELSEIF (ustring(1:4) == 'SCAL' .OR. ustring(1:4) == 'FACT' .OR.  &  
          (informt >= 3 .AND. ustring(1:4) == 'RXFT')) THEN     
     ! ...        error if any one of these formats is used
     IF (kodparm > 0) THEN
        WRITE (fustdout,9005) chrarry(kodparm), fmt(1:7)
        WRITE (fuclog,9005) chrarry(kodparm), fmt(1:7)
9005    FORMAT (/'*** INPUT ERROR: defining ',a,' ***' /tr10,a,  &
             ' input is not valid for this parameter')
        ierr(63) = .TRUE.
     ELSE
        WRITE (fustdout,9015) fmt(1:7)
        WRITE (fuclog,9015) fmt(1:7)
        ierr(64) = .TRUE.
     END IF
  ELSE      ! ... read array with SPECIFIED FORMAT
     kodd = 1
     READ(fuins,fmt) (dumarray(ic),ic=1,nxyz)
     ! ... convert from user specified units to internal units (cgs)
     DO ic = 1,nxyz
        dumarray(ic) = dumarray(ic)*unitfac(kodparm)
     END DO
  END IF
  ! ... Load zero into all cells outside the simulation region
  DO  ic = 1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF(ibc(ik,j) == -1) dumarray(ic) = 0._kdp
  END DO
  DEALLOCATE(dumnx, dumnz)
  RETURN
  ! ...   ---------------- format and error statements -------------------------
200 WRITE (fustdout,9010) chrarry(kodparm)
  WRITE (fustdout,'(tr10,a)') 'Problem with '//fmt(1:7)//' input'
  WRITE (fustdout,9025) nrxtype
  WRITE (fuclog,9010) chrarry(kodparm)
9010 FORMAT (/'***  ERROR defining ', a4, ' ***')
  WRITE (fuclog,'(tr10,a)') 'Problem with '//fmt(1:7)//' input'
  WRITE (fuclog,9025) nrxtype
9025 FORMAT (10X, 'Check that',i2,' rock types are defined')
  ierr(65) = .TRUE.
  DEALLOCATE(dumnx, dumnz)
END SUBROUTINE rd
