SUBROUTINE sink(nppnx,flag)
  ! ... Purpose:  To read source/sink data, and simulation time period data,
  ! ...           and any parameter groups that change
  USE machine_constants, ONLY: kdp
  USE f_units
  USE bc
  USE parameters          !!$*** for flint
  USE control
  USE fdeq, ONLY: ioptupst
  USE i_c
  USE mesh
!!$  USE parameters
  USE source
  USE variables
  USE units
  IMPLICIT NONE
  INTEGER, INTENT(IN OUT) :: nppnx     ! ... number of times SINK has been called
  INTEGER, INTENT(IN OUT) :: flag
  ! 
!!$  CHARACTER(LEN=10), EXTERNAL :: uppercase
  CHARACTER(LEN=4) :: csrcetype, chrtime
  CHARACTER(LEN=76) :: dumstring, duminput, csrceunit1, csrceunit2
  INTEGER :: ic, ichgpr, icxf, icxl, iduma, j, k, kdum, kodd, ldum,  &
       lll, nparms
  INTEGER :: a_err
  LOGICAL :: lchgpr, errflag
  REAL(KIND=kdp) :: delttemp
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: aprint
  INTEGER :: alloc_stat
  INTERFACE
     FUNCTION uppercase(string) RESULT(outstring)
       CHARACTER(LEN=*), INTENT(IN) :: string
       CHARACTER(LEN=LEN(string)) :: outstring
     END FUNCTION uppercase
  END INTERFACE     !!$*** for flint
  INTERFACE         !!$*** for flint
     SUBROUTINE phreg(p,h)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
     END SUBROUTINE phreg
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
  END INTERFACE
  !
  ! ...   ichgpr 1:read in all new print/plot options, 0:no change
  ! ...   lchgpr (logical) T:read in new print/plot options, F:no change
  ! ...   nparams - number of parameters for which new values will be input
  ! ...   nsrce - number of sources; -1:no change (after 1st pumping period)
  ! ...   delt - initial time step
  ! ...         -1: use the time step calculated for the previous step
  ! ...   tchg - time of change in parameters or b.c. values
  ! ...         -1: terminate simulation
  !
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.4 $//$Date: 2007/10/17 18:41:54 $'
  ALLOCATE(aprint(nxyz), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: readarray'
    ierr(179) = .TRUE.
    RETURN
  END IF
  !     ------------------------------------------------------------------
  nppnx = nppnx + 1
  nppn = nppnx
  ! ...     if this is not the first call to sink, set temporary delt
  IF (nppn > 1) THEN
     delttemp = delt
  ELSE
     delttemp = 0._kdp
  END IF
  ! ... read and write simulation period, source, and changed parameter information
  IF (informt == 1) THEN
     READ(fuins,'(1X)')
     READ(fuins,*) tchg, delt, nsrce, lchgpr, nparms
     ichgpr = 0
     IF (lchgpr) ichgpr = 99
  ELSEIF(informt == 2) THEN
     READ(fuins,*,ERR=90) tchg, delt, nsrce, ichgpr, nparms
     IF (ichgpr >= 1) ichgpr = 99
  ELSEIF(informt == 3) THEN
     READ(fuins,*,ERR=90) tchg, delt, nsrce, lchgpr, nparms
     ichgpr = 0
     IF (lchgpr) ichgpr = 1
  END IF
  ! ... Terminate the simulation if TCHG is -1
  IF(tchg <= 0._kdp) THEN
     flag = 99
     DEALLOCATE(aprint)
     RETURN
  END IF
  IF (iyr == 2) THEN
     ! ...  convert time of change from years to seconds (cgs units)
     tchg = sec(tchg)
     IF (delt /= -1.0_kdp) delt = sec(delt)
  END IF
  ! ...         check for input errors
  errflag = .FALSE.
  IF (tchg <= time) THEN
     WRITE(fustdout,9030) '*** ERROR: Simulation period input ***',  &
          'Ending time of simulation period is less than starting time'
     WRITE(fuclog,9030) '*** ERROR: Simulation period input ***',  &
          'Ending time of simulation period is less than starting time'
     ierr(39) = .TRUE.
     errflag = .TRUE.
  END IF
  IF (delt == 0._kdp) THEN     ! ... Must use equal test since -1 is valid
     WRITE(fustdout,9030) '*** ERROR: Simulation period input ***',  &
          'DELT cannot equal 0'
     WRITE(fuclog,9030) '*** ERROR: Simulation period input ***',  &
          'DELT cannot equal 0'
     ierr(38) = .TRUE.
     errflag = .TRUE.
  END IF
  IF (delt > -1.0_kdp .AND. delt < mindelt) THEN 
     WRITE(fustdout,9030) '*** ERROR: Simulation period input ***',  &
          'DELT is less than minimum specified value'
     WRITE(fuclog,9030) '*** ERROR: Simulation period input ***',  &
          'DELT is less than minimum specified value'
     ierr(91) = .TRUE.
     errflag = .TRUE.
  END IF
  IF (delt <= -1._kdp .AND. nppn == 1) THEN
     WRITE(fustdout,9030) '*** ERROR: Simulation period input ***',  &
          'DELT cannot be -1 for first simulation period'
     WRITE(fuclog,9030) '*** ERROR: Simulation period input ***',  &
          'DELT cannot be -1 for first simulation period'
     ierr(41) = .TRUE.
     errflag = .TRUE.
  END IF
  IF(errflag) THEN
    DEALLOCATE(aprint)
    RETURN
  END IF
  ! ...    Use the time step calculated from the previous time
  ! ...       step if DELT=-1 and this is not the first call to SINK
  IF (delt <= -1._kdp .AND. nppn > 1) delt = delttemp
  IF (ioptpr(3) >= 1) WRITE(fustdout,'(/a,i3,a/)') '--- Simulation Period Number', nppn,' ---'
  IF (iyr == 2) THEN
     WRITE(fuclog,9005) '--- Simulation Period Number',nppn,' ---',  &
          'Period begins     ',year(time),' (yr)',  &
          'Period ends       ',year(tchg),' (yr)',  &
          'Initial time step ',year(delt),' (yr)',  &
          'Number of parameters with new data input: ',nparms
     WRITE(fupdef,9005) '--- Simulation Period Number',nppn,' ---',  &
          'Period begins     ',year(time),' (yr)',  &
          'Period ends       ',year(tchg),' (yr)',  &
          'Initial time step ',year(delt),' (yr)',  &
          'Number of parameters with new data input: ',nparms
9005 FORMAT (/a,I3,a/tr10,a,1pg14.7,a/tr10,a,1pg14.7,a/tr10,a,1pg14.7,a/tr10,a,i2)
  ELSE
     WRITE(fuclog,9005) '--- Simulation Period Number',nppn,' ---',  &
          'Period begins     ',time,' (s)',  &
          'Period ends       ',tchg,' (s)',  &
          'Initial time step ',delt,' (s)',  &
          'Number of parameters with new data input: ',nparms
     WRITE(fupdef,9005) '--- Simulation Period Number',nppn,' ---',  &
          'Period begins     ',time,' (s)',  &
          'Period ends       ',tchg,' (s)',  &
          'Initial time step ',delt,' (s)',  &
          'Number of parameters with new data input: ',nparms
  END IF
  IF(informt <= 2) THEN
     IF (ichgpr == 0) THEN
        WRITE(fuclog,'(TR10,A)') 'Print/plot options are not changed'
     ELSE IF(ichgpr > 0) THEN
        WRITE(fuclog,'(TR10,A)') 'New print/plot options are input'
     ELSE
        WRITE(fuclog,'(TR10,I1,A)') ABS(ichgpr), ' print/plot options are changed'
        IF(ABS(ichgpr) > 8) THEN
           WRITE(fustdout,9030) '*** WARNING - Resetting Print/Plot value ***',  &
                'Max number of print options to change limited to 7'
           WRITE(fuclog,9030) '*** WARNING - Resetting Print/Plot value ***',  &
                'Max number of print options to change limited to 7'
           ichgpr = -7
           WRITE(fuclog,'(TR10,I1,A)') ABS(ichgpr), ' print/plot options changed'
        END IF
     END IF
  ELSEIF(informt == 3) THEN
     IF (ichgpr == 0) THEN
        WRITE(fuclog,'(TR10,A)') 'Print/plot options are not changed'
     ELSEIF (ichgpr > 0) THEN
        WRITE(fuclog,'(TR10,A)') 'New print/plot options are input'
     ENDIF
  ENDIF
  !---------------- Sources and Sinks -----------------
  ! ... Allocate space if this is the first call
  IF(nppn == 1) THEN
     CALL alloc_source(a_err)
     IF (a_err /= 0) THEN  
        PRINT *, "Array allocation failed: sink, alloc_source"  
        ierr(199) = .TRUE.
        DEALLOCATE(aprint)
        RETURN
     ENDIF
  ENDIF
  IF (nsrce == -1 .AND. nppn > 1) THEN
     WRITE(fuclog,'(TR10,A)') 'No change in sources or sinks from the previous simulation period'
     WRITE(fupdef,'(TR10,A)') 'No change in sources or sinks from the previous simulation period'
  ELSEIF (nsrce == 0) THEN
     WRITE(fuclog,'(TR10,A)') 'No sources or sinks active during this simulation period'
     WRITE(fupdef,'(TR10,A)') 'No sources or sinks active during this simulation period'
     ! ... Shut off all sources (wells and nodes)
     iq = 0
     q = 0._kdp
     qh = 0._kdp
     qt = 0._kdp
     qhi = 0._kdp
  ELSEIF (nsrce >= 1) THEN          ! ... Active source/sink wells this simulation period
     IF(nsrce == 1) THEN
        WRITE(fuclog,'(TR10,2A,I3)') 'New values input for well or point sources or sinks'
        WRITE(fupdef,'(TR10,2A,I3)') 'New values input for well or point sources or sinks'
     ELSEIF(nsrce == 2) THEN
        WRITE(fuclog,'(TR10,2A,I3)') 'New values input for well and point sources or sinks'
        WRITE(fupdef,'(TR10,2A,I3)') 'New values input for well and point sources or sinks'
     END IF
     ! ... do not re-initialize arrays. unchanged values carry over to new simulation period
     ! ...      read in data for each source/sink that changes this period
     DO  ldum=1,nsrce
        IF (informt == 1) THEN     ! ... original input format, only well input accepted
           CALL sinkwell
           IF(flag > 0) THEN
             DEALLOCATE(aprint)
             RETURN
           END IF
        ELSE
           ! ... read in character string containing input-type and units
           ! ...     and strip leading blanks
           READ(fuins,'(A)') duminput
           duminput = ADJUSTL(duminput)
           ! ... read source/sink input-type and string containing dimension units
           READ(duminput,'(A4,A)') csrcetype, dumstring
           ! ... split unit string into separate components
           iduma = INDEX(dumstring,')')
           csrceunit1 = dumstring(1:iduma)
           csrceunit2 = dumstring(iduma+1:76)   
           csrcetype = uppercase(csrcetype)
           IF (csrcetype == 'WELL') THEN              ! ... well input for source/sinks
              WRITE(fupdef,'(/TR10,A,I2,A)') '--- Well Sources/Sinks ---'
              ! ...            get dimension units for source/sink well input
              CALL getunits('SRC1',csrceunit1,kodd)
              CALL getunits('SRC2',csrceunit2,kodd)
              CALL sinkwell
           ELSE IF (csrcetype == 'NODE') THEN          ! ... node input for source/sinks
              WRITE(fupdef,'(/TR10,A,I2,A)') '--- Point Sources/Sinks ---'
              ! ...            get dimension units for source/sink node input
              CALL getunits('SRC1',csrceunit1,kodd)
              CALL getunits('SRC2',csrceunit2,kodd)
              CALL sinknode
           ELSE                             ! ... error in input type
              WRITE(fustdout,'(/tr10,A/tr15,A/tr15,A)')  &
                   '*** ERROR in source/sink input-type ***',  &
                   'Input-type '//csrcetype//' not recognized.',  &
                   'Acceptable types:  WELL, NODE'
              WRITE(fuclog,'(/tr10,A/tr15,A/tr15,A)')  &
                   '*** ERROR in source/sink input-type ***',  &
                   'Input-type '//csrcetype//' not recognized.',  &
                   'Acceptable types:  WELL, NODE'
              ierr(42) = .TRUE.
              DEALLOCATE(aprint)
              RETURN
           END IF
        END IF
     END DO
     ! ...      display each region slice with source/sink nodes flagged as 1
     WRITE(fupdef,'(/A/)')  '--- Currently Active Source/Sink Cells ---'
     aprint = iq
     CALL printar(2,aprint,ibc,fupdef,1._kdp,10,000)
  END IF
  IF (ichgpr /= 0) THEN
     CALL printopt(ichgpr)        ! ... read in new print/plot options 
  ELSE
     ! ... Reset the print stop flags for unprinted variables to end of time period
     IF(print_press <= 0._kdp) THEN
        timprp=tchg
     ENDIF
     IF(print_enth <= 0._kdp) THEN
        timprh=tchg
     ENDIF
     IF(print_temp <= 0._kdp) THEN
        timprt=tchg
     ENDIF
     IF(print_satn <= 0._kdp) THEN
        timprsat=tchg
     ENDIF
     IF(print_dens <= 0._kdp) THEN
        timprd=tchg
     ENDIF
     IF(print_vis <= 0._kdp) THEN
        timprv=tchg
     ENDIF
     IF(print_pot <= 0._kdp) THEN
        timprpot=tchg
     ENDIF
     IF(print_vel <= 0._kdp) THEN
        timprvel=tchg
     ENDIF
     IF(print_bcflow <= 0._kdp) THEN
        timprbcf=tchg
     ENDIF
     IF(print_source <= 0._kdp) THEN
        timprsrc=tchg
     ENDIF
     IF(print_pmprop <= 0._kdp) THEN
        timprprop=tchg
     ENDIF
     IF(print_poros <= 0._kdp) THEN
        timprpor=tchg
     ENDIF
     IF(print_perm <= 0._kdp) THEN
        timprperm=tchg
     ENDIF
     IF(print_balance <= 0._kdp) THEN
        timprbal=tchg
     ENDIF
     IF(print_dump <= 0._kdp) THEN
        timprdump=tchg
     ENDIF
     IF(print_resid <= 0._kdp) THEN
        timprresid=tchg
     ENDIF
     IF(print_dimno <= 0._kdp) THEN
        timprdimno=tchg
     ENDIF
     IF(print_plotscalar <= 0._kdp) THEN
        timprpscal=tchg
     ENDIF
     IF(print_plotvector <= 0._kdp) THEN
        timprpvec=tchg
     ENDIF
     IF(print_temporal <= 0._kdp) THEN
        timprtemp=tchg
     ENDIF
     ! ... Set the time for next printout
     prtimenxt = MIN(tchg, timprp, timprh, timprt, timprsat, timprd, timprv, timprpot,  &
          timprvel, timprbcf, timprsrc, timprprop, timprpor, timprperm, timprbal, timprdump,  &
          timprresid, timprdimno, timprpscal, timprpvec, timprtemp)
  END IF
  IF (nparms /= 0) THEN           ! ... read new input parameter data blocks
     ! ...        change all KOD's to negative values in order to tell
     ! ...        which parameter arrays have had new values input
     DO  kdum = 1, mptlparm
        kod(kdum) = -kod(kdum)
     END DO
     WRITE(fupdef,'(tr1)')
     CALL readarry(nparms)
     ! ...               check for errors and write new input paramaters
     IF(kod(12) > 0) THEN               ! ...      X-direction node spacing
        WRITE(fuclog,'(A)') 'ERROR - X dimensions of cells cannot be redefined'
        ierr(43) = .TRUE.
     END IF
     IF(kod(13) > 0) THEN               ! ...      Y-direction node spacing
        WRITE(fuclog,'(A)') 'ERROR - Y dimensions of cells cannot be redefined'
        ierr(43) = .TRUE.
        DEALLOCATE(aprint)
        RETURN
     END IF
     IF(kod(14) > 0) THEN               ! ...      Z-direction node spacing
        WRITE(fuclog,'(A)') 'ERROR - Z dimensions of cells cannot be redefined'
        ierr(43) = .TRUE.
        DEALLOCATE(aprint)
        RETURN
     END IF
     ! ... input data for Enthalpy/Temperature and Pressure
     IF (kod(10) > 0 .OR. kod(11) > 0) CALL inputph(.FALSE.)
     ! ... input data for Basal Heat Flux (Conductive)
     IF (kod(9) > 0 .OR. kod(18) > 0) CALL init_wr_bc(.FALSE.)
     ! ... Printout new parameters for rock types
     ! ...         Porosity, Permeability, Thermal Conductivity
     ! ...         Specific Heat, Rock Density, Compressibility
     ! ...     also calc new rock parameters if any are f(Temperature) or f(Time)
     IF (kod(1) > 0 .OR. kod(2) > 0 .OR. kod(3) > 0 .OR.  &
          kod(4) > 0 .OR. kod(5) > 0 .OR. kod(6) > 0 .OR.  &
          kod(7) > 0 .OR. kod(8) > 0 .OR. lrxftreq .OR. lrxftimreq)  &
          CALL init_wr_pmp(.FALSE.,flag)
     ! ... Recalculate LITHOSTATIC PRESSURE if new rock density or new
     ! ...      fluid pressures were input
     IF (kod(7) > 0 .OR. kod(10) > 0) THEN
        CALL preslith
        WRITE(fupdef,9020) 'New Lithostatic Pressure Limit', unitlabel(11)
        CALL printar(2,plith,ibc,fupdef,unitfac(11),24,000)
     END IF
     ! ... Recalculate media and fluid properties using new parameters
     ! ***??? are rockparm call and permftime call needed here???
     CALL tcalc             ! ... Calculate transmissivity and thermal conductivity for the cell faces
     CALL phreg(p,h)        !***this will miss the region for associated enthalpy data
     CALL enthrock          ! ... calculate enthalpy of the matrix
     CALL prpty(0)          ! ... calculate fluid property coefficients and derivatives
     CALL wellallo          ! ... determine the heat and mass allocations of source/sink terms
     IF (ioptupst) CALL upstre     ! ... determine the upstream node for each cell face
     CALL storativ          ! ... calculate mass and energy storativities
     ! ... save mass & energy storativities for restart in case time step is cut
     DO  lll = 1, npiccons
        ic = npic(lll)
        xmoldt(ic) = xm(ic)
        enoldt(ic) = en(ic)
        xmoldt2(ic) = xm(ic)
        enoldt2(ic) = en(ic)
     END DO
     ! ... write message to output files if new values for P,T,or H were input
     IF (kod(10) > 0 .OR. kod(11) > 0) THEN
        IF (kod(11) > 0) THEN
           WRITE(fuclog,9025) '*****  New values input at this time for P *****'
           WRITE(fupdef,9025) '*****  New values input at this time for P *****'
           IF (ioptprta(1) >= 1) &
                WRITE(fup,9025) '*****  New values input at this time for P *****'
9025       FORMAT(//tr5,a/)
        END IF
        IF (kod(10) > 0) THEN
           WRITE(fuclog,9025) '*****  New values input at this time for H or T *****'
           WRITE(fupdef,9025) '*****  New values input at this time for H or T *****'
           IF(ioptprta(2) >= 1)  &
                WRITE(fuen,9025)  '*****  New values input at this time for H or T *****'
           IF(ioptprta(3) >= 1) &
                WRITE(fut,9025) '*****  New values input at this time for H or T *****'
           IF(ioptprta(4) >= 1) &
                WRITE(fusat,9025) '*****  New values input at this time for H or T *****'
        END IF
        IF (ioptpr(6) > 0) THEN      ! ... write new header with dimensional units to temporal plot file
           WRITE(fuplttimser,*)
           IF (iyr == 2) THEN
              chrtime = '(yr)'
           ELSE
              chrtime = '(s) '
           END IF
           ! ...        write units header 
           ! ...           for wide (132 column) output
           IF (ioptpr(2) >= 1) THEN
              WRITE (fuplttimser,'(A,A,A,A)') ' Step   Time   ', '   I  J  K   ',  &
                   'Pressure    Enthalpy  Temperature  Saturation',  &
                   'Water Mass Frac'
              WRITE (fuplttimser,'(8X,A,14X,3A)') chrtime,  &
                   unitlabel(11), unitlabel(10), unitlabel(15)
           ELSE
              ! ...            for narrow (80 column) output
              WRITE(fuplttimser,'(A,A,A)') ' Step   Time   ', '   I  J  K   ',  &
                   'Pressure    Enthalpy  Temperature  Saturation'
              WRITE(fuplttimser,'(8X,A,14X,3A)') chrtime, unitlabel(11),  &
                   unitlabel(10), unitlabel(15)
           END IF
        ENDIF
        IF (ioptpr(5) == 1 .AND. ltime >= 1) WRITE(fusrc,9025)
        IF (ioptpr(5) == 1 .AND. ltime >= 1) WRITE(fubcf,9025)
        ! ... Calculate the water and steam velocity for printing if necessary
        IF (ABS(print_vel) > 0._kdp .OR. ABS(print_plotvector) > 0._kdp .OR.  &
             ABS(print_plotscalar) > 0._kdp .OR.  &
             ioptprta(9) >= 1 .OR. ioptprta(10) >= 1 .OR.  &
             ioptprta(11) >= 1 .OR. ioptprta(15) >= 1 .OR.  &
             ioptprta(16) >= 1 .OR. ioptprta(17) >= 1 .OR. ioptprta(27) >=1) THEN
           CALL velocity
        END IF
        IF (ioptpr(1) == 1) THEN
           CALL plotexpl                ! ... write Explorer plotfile
        ELSEIF(ioptpr(1) == 8 .OR. ioptpr(1) == 9) THEN
           CALL plotb2                  ! ... write B2 plotfile 
        END IF
        IF(ltime > 0) CALL pdata(flag)        ! ...  print results
     END IF
     ! ... change all KOD back to positive values for printing
     DO  kdum = 1, mptlparm
        kod(kdum) = ABS(kod(kdum))
     END DO
  END IF
  WRITE(fupdef,*)
  DEALLOCATE(aprint)
  RETURN
  ! ... error
90 WRITE(fustdout,9015)
  WRITE(fuclog,9015)
9015 FORMAT (/'***** STOP - Error in Simulation Time Period Input *****'/tr10,  &
       'Check that input data agree with variable types')
  ierr(103) = .TRUE.
  flag = flag+1 
  DEALLOCATE(aprint)
  RETURN
9020 FORMAT (/, '--- ', a, ' ---  ', a)
9030 FORMAT (/TR10,a/tr15,a)
END SUBROUTINE sink
