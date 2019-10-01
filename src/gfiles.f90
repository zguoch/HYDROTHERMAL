SUBROUTINE gfiles(nplace,ufile,flag)
  ! ... Purpose: Opens the files for reading input and writing results.
  ! ... Units used:
  ! ...    7 - (input & output) stripped input file (from fuinc) after removing comment lines
  ! ...    8 - (input) Data_file; Main input data file
  ! ...   10 - (output) Out_Probdefine; Information defining the simulation problem
  ! ...   11 - (output) Out_pressure; Pressure fields
  ! ...   12 - (output) Out_enthalpy; Enthalpy fields
  ! ...   13 - (output) Out_temperature; Temperature fields
  ! ...   14 - (output) Out_saturation; Saturation fields
  ! ...   15 - (output) Out_density; Fluid density fields
  ! ...   16 - (output) Out_viscosity; Fluid viscosity fields
  ! ...   17 - (output) Out_potential; Potential and potential head fields
  ! ...   18 - (output) Out_velocity; Velocity fields
  ! ...   19 - (output) Out_bcflow; Boundary condition flow rates
  ! ...   20 - (output) Out_source; Source flow rates
  ! ...   21 - (output) Out_pmthermalprop; Porous media physical and thermal properties
  ! ...   22 - (output) Out_porosity; Porosity field
  ! ...   23 - (output) Out_permeability; Permeability and relative permeability fields
  ! ...   24 - (output) Out_balance; Global balance tables for flow and heat
  ! ...   25 - (output) Out_restartdump; Dump of parameters for simulation restart
  ! ...   26 - (output) Out_residual; Mass and energy residual fields
  ! 
  ! ...   27 - (output) Plot_scalar; Plotfile for scalar data
  ! ...   28 - (output) Plot_vector; Plotfile for vector data
  ! ...   29 - (output) Plot_timeseries; Plotfile for temporal series data
  !
  ! ...   30 - (output) C_log; calculation log with time step and iteration data
  !
  ! ...   32 - (output) Out_dimensionless; Dimensionless number fields
  ! ...   35 - (output) Plot_scalar2; Plotfile number 2 for scalar data
  USE f_units
  USE bc, ONLY: unconfined
  USE control
  USE units
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nplace       ! ...  location to go in subroutine to open appropriate files
  CHARACTER(LEN=80), INTENT(IN OUT) :: ufile     ! ... filename
  INTEGER, INTENT(OUT) :: flag        ! ... status flag
  !
  CHARACTER(LEN=80) :: upltdir
  CHARACTER(LEN=5) :: chrdum1
  CHARACTER(LEN=11) :: fmt
  CHARACTER(LEN=31) :: name
  INTEGER :: ios, length
  LOGICAL :: lex, lerror
  INTERFACE
     FUNCTION uppercase(string) RESULT(outstring)
       CHARACTER(LEN=*), INTENT(IN) :: string
       CHARACTER(LEN=LEN(string)) :: outstring
     END FUNCTION uppercase
  END INTERFACE
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.29 $//$Date: 2007/12/20 20:32:47 $'
  !     ------------------------------------------------------------------
  !...
  flag = 0 
  lerror = .FALSE.
  !!$ luf = lnblnk(ufile)
  IF (nplace == 1) THEN
     ! ... open main input (fuins) and output (fupdef) units file
     ! ... open unit (fuinc) for commented data file
     upltdir = ' '
     IF(.NOT.gui) THEN
        ! ... For HT standalone execution
11      WRITE(fustdout,*) 'Enter name of input data file'
        READ(fustdin,'(A)') name
        name = ADJUSTL(name)
        IF(uppercase(name) == 'Q') THEN
           STOP 'Simulation aborted'
        END IF
        fname = name
        INQUIRE(FILE=fname,EXIST=lex)
        IF(lex) THEN
           OPEN(fuinc,FILE=fname,IOSTAT=ios,STATUS='OLD',POSITION='REWIND',ACTION='READ')
           IF (ios > 0) THEN    
              lerror = .TRUE.
              WRITE(*,*) 'ERROR: Error opening file ', fname
           ENDIF
           ufile=fname
        ELSE
           WRITE(fustdout,'(2A)') TRIM(fname),' DOES NOT EXIST'
           GO TO 11
        END IF
        WRITE(fustdout,*) 'Enter i.d. for output files (up to 14 characters)'
        READ(fustdin,'(A)') name
        name = ADJUSTL(name)
        IF(uppercase(name) == 'Q') THEN
           STOP 'Simulation aborted'
        END IF
        IF(name(15:15) /= ' ') WRITE(fustdout,*)  &
             'Filenames will be truncated to 20 characters'
        usuff = '.'//TRIM(name(1:14))
     ELSEIF(gui) THEN
        ! ... For HT execution under GUI
        ! ...   input filename defined by GUI as 'ht.in'
        INQUIRE(FILE=fname,EXIST=lex)
        IF (lex) THEN
           OPEN(fuinc,FILE=fname,IOSTAT=ios,STATUS='OLD',POSITION='REWIND',ACTION='READ')
          IF (ios > 0) THEN
              WRITE(*,*) 'ERROR: Error opening file ', fname
              flag = 5
              RETURN
           ENDIF
           ufile = fname
        ELSE
           ierr(2) = .TRUE.
           RETURN
        END IF
        ! ... usuff is set by GUI as ' '
        usuff = ' '
        fustdout = 66              ! ... Makes standard out another disc file
        OPEN(fustdout,FILE='Stdout.HT',IOSTAT=ios,STATUS='REPLACE',ACTION='WRITE')
        IF (ios > 0) THEN
           lerror = .TRUE.
           WRITE(*,*) 'ERROR: Error opening file ', fname
           flag = 6
           RETURN
        ENDIF
        ! ... Alternative would be to redirect to dev/null on UNIX
     END IF
     !$$  luf = lnblnk(ufile)
     ! ... unit for input file after stripping out comment statements
     fname = 'Dat.strip'
     OPEN(fuins,FILE=fname,IOSTAT=ios,STATUS='REPLACE')
     IF (ios > 0) lerror = .TRUE.
     ! ... Open general problem definition output file
     fname='Out_Probdefine'//usuff
     OPEN(fupdef,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     ! ... Open calculation log file
     fname='Calc_log'//usuff
     OPEN(fuclog,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     IF(lerror) THEN
        WRITE(*,*) 'ERROR: Error opening one or more output files in gfiles'
        flag = 6
        RETURN
     END IF
  ELSE IF (nplace == 2) THEN
     ! ... open IC_pressporo for input of initial pressure and porosity
     fname = TRIM(ufile)
     OPEN(fuicpp, FILE=fname, STATUS='old', IOSTAT=ios)
     IF (ios > 0) THEN
        lerror = .TRUE.
        WRITE(*,*) 'ERROR: Error opening file ', fname
        flag = 5
        ierr(3) = .TRUE.
        RETURN
     END IF
  ELSE IF (nplace == 3) THEN
     ! ... Open print output files
     fname='Out_pressure'//usuff
     OPEN(fup,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_enthalpy'//usuff
     OPEN(fuen,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_temperature'//usuff
     OPEN(fut,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_saturation'//usuff
     OPEN(fusat,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_watertable'//usuff
     OPEN(fuwtelev,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_density'//usuff
     OPEN(fuden,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_viscosity'//usuff
     OPEN(fuvis,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_potential'//usuff
     OPEN(fupot,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_velocity'//usuff
     OPEN(fuvel,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_bcflow'//usuff
     OPEN(fubcf,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_bcflow2'//usuff
     OPEN(fubcf2,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_source'//usuff
     OPEN(fusrc,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_pmthermprop'//usuff
     OPEN(futhp,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_porosity'//usuff
     OPEN(fupor,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_permeability'//usuff
     OPEN(fuperm,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_balance'//usuff
     OPEN(fubal,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
!!$     fname='Out_restartdump'//usuff
!!$     OPEN(fuorst,FILE=fname,IOSTAT=ios,ACTION='WRITE')
!!$     IF (ios > 0) lerror = .TRUE.
     fname='Out_residual'//usuff
     OPEN(fures,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     fname='Out_dimensionless'//usuff
     OPEN(fudimno,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     IF(lerror) THEN
        WRITE(*,*) 'ERROR: Error opening one or more output files in gfiles '
        flag = 6
        RETURN
     END IF
     IF (iyr == 2) THEN
        chrdum1 = '(yr)'
     ELSE
        chrdum1 = '(s)'
     END IF
     ! ...  write header to file for output flux data for source and constant value nodes
     IF (ioptpr(2) == 0) THEN               ! ...        80 columns
        WRITE(fusrc, '(38X,A)') '(+ = flux in;  - = flux out of region)'
        WRITE(fusrc, '(A,A)') 'Step   Time     I  J   K  Type   ',  &
             'MassFlo TotHeatFlo  AdvcFlo   CondFlo   Advec H'
        WRITE(fusrc, '(A,A,A,A,A)') '       ', chrdum1,  &
             '                     ', ' (g/s)    (erg/s)   (erg/s)   (erg/s)   (erg/g)'
     ELSE                                   ! ...         132 or 159 columns
        WRITE(fusrc, '(44X,A)')  &
             '(+ = flux in;  - = flux out of region)   (# adj nodes)'
        WRITE(fusrc, '(A,A)') ' Step     Time     I   J   K   Type      '  &
             , 'MassFlow  Tot.HeatFlo   AdvecFlo    CondFlo     Advec H'
        WRITE(fusrc, '(A,A,A,A,A)') '         ', chrdum1,  &
             '                           ',  &
             '  (g/s)     (erg/s)     (erg/s)     (erg/s)     (erg/g)'
     END IF
     ! ...      write header to file for output flux data for boundary cells
     IF (ioptpr(2) == 0) THEN               ! ...        80 columns
        WRITE(fubcf, '(38X,A)') '(+ = flux in;  - = flux out of region)'
        WRITE(fubcf, '(A,A)') 'Step   Time     I  J   K  Type   ',  &
             'MassFlo TotHeatFlo  AdvcFlo   CondFlo   Advec H'
        WRITE(fubcf, '(A,A,A,A,A)') '       ', chrdum1,  &
             '                     ', ' (g/s)    (erg/s)   (erg/s)   (erg/s)   (erg/g)'
     ELSE                                  ! ...         132 or 159 columns
        WRITE(fubcf, '(44X,A)')  &
             '(+ = flux in;  - = flux out of region)   (# adj nodes)'
        WRITE(fubcf, '(A,A)') ' Step     Time     I   J   K   Type      '  &
             , 'MassFlow  Tot.HeatFlo   AdvecFlo    CondFlo     Advec H'
        WRITE(fubcf, '(A,A,A,A,A)') '         ', chrdum1,  &
             '                           ',  &
             '  (g/s)     (erg/s)     (erg/s)     (erg/s)     (erg/g)'
     END IF
  ELSE IF (nplace == 4) THEN
     ! ... open file fuicpp for output of initial pressure and porosity data
     fname='IC_pressporo'//usuff
     OPEN(fuicpp, FILE=fname, IOSTAT=ios, ACTION='WRITE')
     IF (ios > 0) THEN
        lerror = .TRUE.
        WRITE(*,*) 'ERROR: Error opening file ', fname
        flag = 6
        RETURN
     ENDIF
  ELSE IF (nplace == 5) THEN
     ! ... Open fuorst for output of current P and H and T data for restart file
     fname='Out_restartdump'//usuff
     OPEN(fuorst, FILE=fname, IOSTAT=ios, ACTION='WRITE')
     IF (ios > 0) THEN
        lerror = .TRUE.
        WRITE(*,*) 'ERROR: Error opening file ', fname
        flag = 6
        RETURN
     ENDIF
  ELSE IF (nplace == 6) THEN
     !  ------ open for output of explorer plotfiles -------
     fname='Plot_scalar'//usuff
     OPEN(fupltsca,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     IF(unconfined) THEN
        fname='Plot_scalar2'//usuff
        OPEN(fupltsca2,FILE=fname,IOSTAT=ios,ACTION='WRITE')
        IF (ios > 0) lerror = .TRUE.
     END IF
     fname='Plot_vector'//usuff
     OPEN(fupltvec,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     IF(lerror) THEN
        WRITE(*,*) 'ERROR: Error opening one or more output files in gfiles '
        flag = 6
        RETURN
     END IF
  ELSE IF (nplace == 7) THEN
     ! ... open  for output of gnuplot or IDL or columnar xyz plotfile
     ! ...      open formatted or unformatted plotfile for gnuplot or IDL
     SELECT CASE (ioptpr(1))
     CASE (2)         !... open unformatted plotfile for gnuplot
        fmt='unformatted'
     CASE (3)         !... open formatted plotfile for gnuplot
        fmt='formatted'
     CASE (4)         !... open unformatted plotfile for IDL
        fmt='unformatted'
     CASE (5)         !... open formatted plotfile for IDL
        fmt='formatted'
     CASE (6)         !... open formatted plotfile for columnar
        fmt='formatted'
     CASE (8)         !... open formatted plotfile for Basin2
        fmt='formatted'
     CASE (9)         !... open unformatted plotfile for Basin2
        fmt='unformatted'
     END SELECT
     fname='Plot_scalar'//usuff
     OPEN(fupltsca,FILE=fname,FORM=fmt,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     IF(unconfined) THEN
        fname='Plot_scalar2'//usuff
        OPEN(fupltsca2,FILE=fname,FORM=fmt,IOSTAT=ios,ACTION='WRITE')
        IF (ios > 0) lerror = .TRUE.
     END IF
     fname='Plot_vector'//usuff
     OPEN(fupltvec,FILE=fname,FORM=fmt,IOSTAT=ios,ACTION='WRITE')
     IF (ios > 0) lerror = .TRUE.
     IF(lerror) THEN
        WRITE(*,*) 'ERROR: Error opening one or more plot output files in gfiles '
        flag = 6
        RETURN
     END IF
  ELSE IF (nplace == 8) THEN
     !  ------ open  for output of time series to plotfile -------
     fname='Plot_timeseries'//usuff
     OPEN(fuplttimser,FILE=fname,IOSTAT=ios,ACTION='WRITE')
     IF(ios > 0) THEN
        lerror = .TRUE.
        WRITE(*,*) 'ERROR: Error opening file ', fname
        flag = 6
        RETURN
     END IF
     !  ------ write header to time series output file
     IF (iyr == 2) THEN
        chrdum1 = '(yr)'
     ELSE
        chrdum1 = '(s)'
     END IF
     ! ...        write units header 
     IF (ioptpr(2) >= 1) THEN         ! ... for wide (132 column) output
        WRITE (fuplttimser,'(A,A,A,A)') ' Step   Time   ', '   I  J  K   ',  &
             'Pressure    Enthalpy  Temperature  Saturation ',  &
             'Water Mass Frac'
        WRITE (fuplttimser,'(tr8,a,tr14,3a)') chrdum1,  &
             unitlabel(11), unitlabel(10), unitlabel(15)
     ELSE                             ! ... for narrow (80 column) output
        WRITE(fuplttimser,'(A,A,A)') ' Step   Time   ', '   I  J  K   ',  &
             'Pressure    Enthalpy  Temperature  Saturation'
        WRITE(fuplttimser,'(tr8,a,tr14,3a)') chrdum1,  &
             unitlabel(11), unitlabel(10), unitlabel(15)
     END IF
  ELSE        ! ...  invalid value of NPLACE
     WRITE (fustdout,'(tr5,a)') '*** ERROR: Invalid value of NPLACE in gfiles ***'
     WRITE (fuclog,'(tr5,a)') '*** ERROR: Invalid value of NPLACE in gfiles ***'
     ierr(4) = .TRUE.
     RETURN
  END IF
END SUBROUTINE gfiles
