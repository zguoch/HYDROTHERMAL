SUBROUTINE inputph(initial)
  ! ... Set up and print the initial or new enthalpy, temperature
  ! ...     and pressure arrays
  USE f_units
  USE bc
  USE control
  USE mesh, ONLY: nx, nz, nxyz
  USE units
  USE variables
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: initial
  ! ...   INITIAL - .TRUE.= initial condition values input
  ! ...             .FALSE.= new values input, usually b.c. values
  ! ... Initial is for initial conditions for pressure, enthalpy, temperature
  ! ... otherwise it is boundary conditions for pressure, enthalpy, temperature
  !
  INTEGER :: i, ic, ik, j, k
  INTERFACE
     SUBROUTINE printar(ndim,array,lprnt,fu,cnv,jfmt,nnoppr)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: ndim
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: array
       INTEGER, DIMENSION(:,:), INTENT(IN) :: lprnt
       INTEGER, INTENT(IN) :: fu
       REAL(KIND=kdp), INTENT(IN) :: cnv
       INTEGER, INTENT(IN) :: jfmt
       INTEGER, INTENT(IN) :: nnoppr
     END SUBROUTINE printar
  END INTERFACE     !!$*** for flint
  INTERFACE         !!$*** for flint
     SUBROUTINE enthtemp(p,tc,h)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: p, tc, h
     END SUBROUTINE enthtemp
     SUBROUTINE phreg(p,h)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
     END SUBROUTINE phreg
     SUBROUTINE press(ndenw,h,tc,p)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: ndenw
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: h, tc, p
     END SUBROUTINE press
     SUBROUTINE tempdegc(p,h,tc)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
       REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: tc
     END SUBROUTINE tempdegc
     SUBROUTINE presboil(p,h,tc,ii,jj)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: p, h, tc
       INTEGER, INTENT(IN) :: ii, jj
     END SUBROUTINE presboil
     SUBROUTINE wrnew(funit,dumarray,kodd,convfac)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: funit
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: dumarray
       INTEGER, INTENT(IN) :: kodd
       REAL(KIND=kdp), INTENT(IN OUT) :: convfac
     END SUBROUTINE wrnew
  END INTERFACE
  !     ------------------------------------------------------------------
  ! ...   Kod() =1 - array input as variable value array
  ! ...         =2 - array input as constant value array
  ! ...         =3 - array input as values constant for each row
  ! ...         =4 - array input as values constant for each column
  ! ...         =5 - array input as rocktype values
  ! ...         =50 - special code for writing Explorer plotfiles
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
  ! ...   Kodt  =0 if enthalpy data are input for all nodes
  ! ...          4 if enthalpy data are input for only a subset of nodes
  ! ...          1 if temperature data are input for all nodes
  ! ...          2 if temperature data are input for only a few nodes
  ! ...          3 if temperature data are input with boiling point options
  !
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.12 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  !***todo: modify to print out only the bc values for new values not the whole field
  ! ...                        ENTHALPY/TEMPERATURE
  IF (kod(10) > 0) THEN
     IF (initial .AND. kodt == 4) THEN    ! Check for error
        WRITE(fustdout,9021) 
        WRITE(fuclog,9021)
9021    FORMAT (/10X,'*** ERROR: Problem with initial conditions ***',  &
                /15X,'Enthalpy values not specified for all nodes')
        ierr(16) = .TRUE.
     ELSE IF (initial .AND. kodt == 2) THEN    ! Check for error
        WRITE(fustdout,9020)
        WRITE(fuclog,9020)
9020    FORMAT (/10X,'*** ERROR: Problem with initial conditions ***',  &
                /15X,'Temperature values not specified for all nodes')
        ierr(16) = .TRUE.
     ELSE IF (kodt == 0 .OR. kodt == 4) THEN          ! ... enthalpy values input
        IF (kodt == 4) THEN    ! new enthalpy values for a subset of nodes
           ! copy new enthalpy values from ehassoc to h
           DO ic=1,nxyz
              CALL ictoijk(ic,i,j,k,nx,nz)
              ik = (k-1)*nx + i
              IF (prflag(ik,j) == 1) THEN
                 h(ic) = ehassoc(ic)
              END IF
           END DO
        ELSE    ! initial or new enthalpy values for all nodes
           ! copy all elements in ehassoc into h
           h = ehassoc
        END IF
        IF (initial) THEN 
            WRITE(fupdef,9005) 'Initial Enthalpy', unitlabel(10)
9005       FORMAT (/'--- ', a, ' ---  ', a)
        ELSE
            WRITE(fupdef,9005) 'New values for Enthalpy', unitlabel(10)
        END IF
        IF(kod(10) == 1) THEN
          CALL printar(2,h,ibc,fupdef,unitfac(10),24,000)
        ELSE
          CALL wrnew(fupdef,h,kod(10),unitfac(10))
        END IF
        IF (kod(11) == 60) CALL press(2,h,tc,p)    ! ... calculate hydrostatic pressure gradient
        CALL phreg(p,h)
        CALL tempdegc(p,h,tc)    ! ... calculate the temperature
        IF (initial) THEN
           WRITE(fupdef,9005) 'Initial Temperature', unitlabel(15)
           WRITE(fupdef,'(tr10,a)') 'Computed from initial enthalpy and pressure'
        ELSE
           WRITE(fupdef,9005) 'New Temperature', unitlabel(15)
           WRITE(fupdef,'(tr10,a)') 'computed from new enthalpy and pressure'
        END IF
        IF (kod(10) == 77) WRITE(fupdef,'(tr10,a)') 'Some values on saturated water curve'
        CALL printar(2,tc,ibc,fupdef,unitfac(15),12,000)
        tcassoc = tc
     ELSEIF (kodt == 1 .OR. kodt == 2) THEN        ! ... temperature values input
        IF (kodt == 2) THEN    ! new temperature values for a subset of nodes
           ! copy new temperature values from tcassoc to tc
           DO ic=1,nxyz
              CALL ictoijk(ic,i,j,k,nx,nz)
              ik = (k-1)*nx + i
              IF (prflag(ik,j) == 1) THEN
                  tc(ic) = tcassoc(ic)
              END IF
           END DO
        ELSE    ! initial or new temperatures for all nodes
           ! copy all elements from tcassoc to tc
           tc = tcassoc
        END IF
        IF (kodt == 1) THEN
           IF (initial) THEN
              WRITE(fupdef,9005) 'Initial Temperature', unitlabel(15)
           ELSE
              WRITE(fupdef,9005) 'New values for Temperature', unitlabel(15)
           END IF
           CALL printar(2,tc,ibc,fupdef,unitfac(15),12,000)
        END IF
        CALL enthtemp(p,tc,h)   ! ... calculate enthalpy from temperature and pressure 
        kod(10) = 1
        IF (initial) THEN
           WRITE(fupdef,9005) 'Initial Enthalpy', unitlabel(10)
           WRITE(fupdef,'(tr10,a)') 'Computed from initial temperature and pressure'
        ELSE
           WRITE(fupdef,9005) 'New Enthalpy', unitlabel(10)
           IF (kod(10) == 77) THEN
              WRITE(fupdef,'(tr10,a)') 'Some values on saturated water curve'
           ELSE
              WRITE(fupdef,'(tr10,a)') 'Computed from new temperature and pressure'
           END IF
        END IF
        CALL printar(2,h,ibc,fupdef,unitfac(10),24,000)
        IF (kodt == 2) THEN
           ! ... recalculate the temperature for all nodes using the new P&H values
           CALL phreg(p,h)           ! ... determine the phase region for each node
           CALL tempdegc(p,h,tc)        ! ... calculate temperature
           WRITE(fupdef,9005) 'New Temperatures recalculated', unitlabel(15)
           WRITE(fupdef,'(tr10,a)') 'Computed from new enthalpy and pressure'
           IF (kod(10) == 77) WRITE(fupdef,'(tr10,a)') 'Some values on saturated water curve'
           CALL printar(2,tc,ibc,fupdef,unitfac(15),12,000)
        END IF
        ehassoc = h     ! needed for nodes that are "specified pressure with associated temperature or enthalpy"
     ELSEIF (kodt == 3) THEN         ! ...  boiling point with depth curve requested
        IF (kod(10) == 70) THEN
           IF(kod(11) == 60) THEN
              kod(10) = 1
!!$              !...***later mod for call to presboil
!!$              !..            do 10 mt=1,nxy
!!$              !..               ptopx(mt)=-1.d0
!!$              !..               m=(nz-1)*nxy+mt
!!$              !..               call mtoijk(m,i,j,k,nx,ny)
!!$              !..            ic = (nz-1)*Nx + i + (j-1)*Nx*Nz
!!$              !..            ptopx(mt)=P(ic)
!!$              !..   10       continue
              IF (initial) THEN
                 CALL presboil(p,h,tc,0,0)
                 WRITE(fupdef,9005) 'Initial Temperature', unitlabel(15)
                 WRITE(fupdef,'(tr10,a)') 'for Boiling Point with Depth'
                 CALL printar(2,tc,ibc,fupdef,unitfac(15),12,000)
                 WRITE(fupdef,9005) 'Initial Enthalpy', unitlabel(10)
                 WRITE(fupdef,'(tr10,a)') 'for Boiling Point with Depth'
                 CALL printar(2,h,ibc,fupdef,unitfac(10),24,000)
                 ehassoc = h
                 tcassoc = tc
              ELSE
                 CALL presboil(p,ehassoc,tcassoc,0,0)
                 WRITE(fupdef,9005) 'New Temperature', unitlabel(15)
                 WRITE(fupdef,'(tr10,a)') 'for Boiling Point with Depth'
                 CALL printar(2,tcassoc,ibc,fupdef,unitfac(15),12,000)
                 WRITE(fupdef,9005) 'New Enthalpy', unitlabel(10)
                 WRITE(fupdef,'(tr10,a)') 'for Boiling Point with Depth'
                 CALL printar(2,ehassoc,ibc,fupdef,unitfac(10),24,000)
                 h = ehassoc
                 tc = tcassoc
              END IF
           ELSE 
              ! ... error if boiling point with depth curve requested, but
              ! ...     hydrostatic pressure gradient not requested
              WRITE(fustdout,'(/a)') 'ERROR: Hydrostatic pressure gradient must be selected to '//  &
                   'calculate boiling point with depth'
              WRITE(fuclog,'(/a)') 'ERROR: Hydrostatic pressure gradient must be selected to '//  &
                   'calculate boiling point with depth'
              ierr(44) = .TRUE.
           END IF
        ELSEIF (kod(10) == 75 .AND. (.NOT.initial)) THEN
           WRITE(fupdef,9005) 'New Temperature', unitlabel(15)
           WRITE(fupdef,'(tr10,a)') 'for Boiling Point with Depth for some columns'
           CALL printar(2,tcassoc,ibc,fupdef,unitfac(15),12,000)
           WRITE(fupdef,9005) 'New Enthalpy', unitlabel(10)
           WRITE(fupdef,'(tr10,a)') 'for Boiling Point with Depth for some columns'
           CALL printar(2,ehassoc,ibc,fupdef,unitfac(10),24,000)
           h = ehassoc
           tc = tcassoc
        END IF
     END IF
  END IF
  ! ...                            PRESSURE
  IF (kod(11) > 0) THEN
     ! ...        stop program if pressure is redefined for some temperature options
     ! ...            if boiling point with depth is specified for some columns
     ! ...            (TEMP CALC option 33) then pressure should not be redefined
     IF (kod(10) == 75) THEN
        WRITE(fustdout,9030)
        WRITE(fuclog,9030)
9030    FORMAT (/5X,  &
             '*** ERROR: Pressure should not be redefined when TEMPERATURE '//  &
             'CALC option 33 is used ***'/8X,  &
             'Hydrostatic pressure is already specified for some '//  &
             'columns with this option')
        ierr(17) = .TRUE.
     ELSEIF (kod(10) == 77) THEN
        ! ...        boiling point is specified for some nodes
        ! ...        (TEMP CALC option 33) then pressure should not be redefined
        WRITE(fustdout,9035)
        WRITE(fuclog,9035)
9035    FORMAT (/10X, '*** ERROR: Pressure must not be redefined'/  &
             15X,'when assigning some nodes to boiling temperature **')
        ierr(45) = .TRUE.
     END IF
     ! ...            write initial or new input values for Pressure
     IF (kod(11) /= 60) THEN
        IF (kodt /= 3) THEN
           IF (initial) THEN
              WRITE(fupdef,9005) 'Initial Pressure', unitlabel(11)
           ELSE
              WRITE(fupdef,9005) 'New values for Pressure', unitlabel(11)
           END IF
           IF(kod(11) == 1) THEN
              CALL printar(2,p,ibc,fupdef,unitfac(11),24,000)
           ELSE
              CALL wrnew(fupdef, p, kod(11), unitfac(11))
           END IF
        END IF
     ELSE          ! ...  Hydrostatic conditions
        ! ...     write value for Pressure at top of domain
!!$        IF (initial) THEN
!!$           WRITE(fupdef,9005) 'Pressure at top of problem domain', unitlabel(11)
!!$        ELSE
!!$           WRITE(fupdef,9005) 'New Pressure at top of problem domain',  &
!!$                unitlabel(11)
!!$        END IF
!!$        WRITE(fupdef,9050) 'Uniform Value: ',ptop/unitfac(11)
!!$9050    FORMAT (10X,a,1PG11.4)
        ! ... Compute and print hydrostatic pressure if boiling point
        ! ...      with depth is NOT specified
        IF (kodt /= 3) THEN
           IF (initial) THEN
              CALL press(2,h,tc,p)
              WRITE(fupdef,9005) 'Initial Calculated Hydrostatic Pressure',  &
                   unitlabel(11)
           ELSE
              CALL press(2,ehassoc,tcassoc,p)   
              WRITE(fupdef,9005) 'New Calculated Hydrostatic Pressure',  &
                   unitlabel(11)
           END IF
           kod(11) = 1
           CALL printar(2,p,ibc,fupdef,unitfac(11),24,000)
        ELSEIF (kodt == 3) THEN
           ! ... print pressure for boiling point with depth;
           ! ...      hydrostatic pressure gradient requested
           IF (initial) THEN
              WRITE(fupdef,9005) 'Initial Pressure', unitlabel(11)
              WRITE(fupdef,'(tr10,a)') 'for Boiling Point with Depth'
           ELSE
              WRITE(fupdef,9005) 'New Pressure', unitlabel(11)
              WRITE(fupdef,'(tr10,a)') 'for Boiling Point with Depth'
           END IF
           kod(11) = 1
           CALL printar(2,p,ibc,fupdef,unitfac(11),24,000)
        END IF
     END IF
  ELSE IF (kod(10) == 75 .AND. .NOT.initial) THEN
     ! ...         pressure not specifically redefined
     ! ...            if boiling point with depth is specified for some columns
     ! ...            (TEMP CALC option 33) then pressure is redefined automatically
     ! ...            for those nodes
     kod(10) = 1
     IF (kod(10) == 75) THEN
        WRITE(fupdef,9005) 'New Pressure', unitlabel(11)
        WRITE(fupdef,'(tr10,a)') 'for Boiling Point with Depth for some columns'
        !..          CALL WRNEW(fupdef, P, 1, Unitfac(11))
        CALL printar(2,p,ibc,fupdef,unitfac(11),24,000)
     END IF
  ELSE IF (kod(10) >= 70) THEN
     IF (kod(10) == 77) THEN
        kod(10) = 1
     ELSE
        WRITE(fustdout,9040) ' **** Problem with kod(10) and kod(11) in INPUTPH'
        WRITE(fuclog,9040) ' **** Problem with kod(10) and kod(11) in INPUTPH'
9040    FORMAT (/a)
     END IF
  END IF
  IF (kod(10) > 5 .OR. kod(11) > 5) THEN
     WRITE(fustdout,9045) kod(10), kod(11)
     WRITE(fuclog,9045) kod(10), kod(11)
9045 FORMAT (/' **** ERROR: Problem with KOD values in INPUTPH'/  &
          5X,'KOD(10)=', i3, 'KOD(11)=', i3)
  END IF
END SUBROUTINE inputph
