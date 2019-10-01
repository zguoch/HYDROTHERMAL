SUBROUTINE time_step_ht(flag)
  ! ... Takes the simulation through one time step, including time step
  ! ...      reduction if necessary
  ! ... Sends back a non-zero flag if ending time has been reached 
  ! ...      or an error occurs
  USE f_units
  USE bc
  USE parameters          !!$*** for flint
  USE control
  USE fdeq
  USE i_c, ONLY: plith
  USE mesh
!!$  USE parameters
  USE source
  USE units
  USE variables
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: flag
  !
  CHARACTER (LEN=1) :: aconvrg, atmcntl, atmcut
  CHARACTER (LEN=13) :: uconvrg
  LOGICAL :: ioptlitho, lprint
  INTEGER :: i, ic, ich, ichdiffmax, ichdiffpc, icp, icpdiffmax,  &
       icpdiffpc, idum4, idum5, idum6, idum7, idum8, idum9,  &
       iduma, idumb, idumc, ie, ihchg, ihcon
  INTEGER :: ihpc, ipchg, ipcon, ippc,  &
       iwschg, j, jhchg, jhcon, jhpc, jpchg, jpcon, jppc, jwschg,  &
       k, khchg, khcon, khpc
  INTEGER :: kpchg, kpcon, kppc, kwschg, lconvrg, lll, ntimecut
  INTEGER  nswitch, lnrichk1, lnrichk2
  INTEGER ij, ik
  REAL(KIND=kdp) :: awschgocr, apchgocr, ahchgocr, deltoldt2, dum1, dum2, dum3, fac, &
       gveqnm, gveqne, &
       hchgocr,  &
       hconvrgmx, hconvrgpc, lqir, pchgocr, pconvrgmx,  &
       pconvrgpc, qfluxsp, qhfluxsp,  &
       scrfir,  &
       timeold, upri, vpir, wschgocr
  REAL(KIND=kdp), dimension(3), parameter :: hh=(/0.01e9_kdp, 0.5e9_kdp, 52.0e9_kdp/),  &
       pp=(/0.5e6_kdp, 2.4e8_kdp, 1.0e10_kdp/)
  REAL(KIND=kdp) :: hmin, hmax, pmin, pmax 
  INTERFACE
     SUBROUTINE phreg(p,h)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
     END SUBROUTINE phreg
     SUBROUTINE tempdegc(p,h,tc)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
       REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: tc
     END SUBROUTINE tempdegc
  END INTERFACE
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.4 $//$Date: 2007/10/16 21:17:22 $'
  !     ------------------------------------------------------------------
  ! ... Begin a Time Step 
  flag = 0
  nswitch = 0             !     nswitch is not used
  ltime = ltime + 1
  IF(ltime > ltimemax) THEN
     ierr(190) = .TRUE.
     flag = 2
     RETURN
  ENDIF
  IF (iyr == 2) THEN
     WRITE (fustdout,9005) '--- Begin Time Step No. ',ltime,  &
          ' ---  Time Step Length: ',year(delt),' (yr)'
     IF (ioptpr(3) >= 1) THEN
        WRITE (fuclog,9005) '--- Begin Time Step No. ',ltime,  &
             ' ---  Time Step Length: ',year(delt),' (yr)'
9005    FORMAT(/a,i5,a,1pg12.6,a)
     END IF
  ELSE
     WRITE (fustdout,9005) '--- Begin Time Step No. ',ltime,  &
          ' ---  Time Step Length: ',delt,' (s)'
     IF (ioptpr(3) >= 1) THEN
        WRITE (fuclog,9005) '--- Begin Time Step No. ',ltime,  &
             ' ---  Time Step Length: ',delt,' (s)'
     END IF
  END IF
  ! ... Write which criterion determined the size of the current time step
  CALL timestep_criterion(ltmcntl, atmcntl)
  IF(ierr(125)) go to 900     ! ... error in time step criterion  !---this should not happen, unless there is a bug in the program.
  ! ...     initialize variables for time step
  tmcut = .FALSE.
  atmcut = ' '
  ltmcntl = 0
  ntimecut = 0
  ! ... Save old time level
  timeold = time
  ! ... Set new time level
  time = time + delt
  ! ...     save mass & energy storativities & saturation from last time step
!!$  icall = 0
  DO  lll = 1, npiccons
     ic = npic(lll)
     indoldt(ic) = ind(ic)
     swoldt(ic) = satnw(ic)
     poldt(ic) = p(ic)
     holdt(ic) = h(ic)
     xmoldt2(ic) = xmoldt(ic)
     enoldt2(ic) = enoldt(ic)
     xmoldt(ic) = xm(ic)
     enoldt(ic) = en(ic)
  END DO
!!$  ! ...      return to here if time step is cut or seepage cell switched.
  ! ...      return to here if time step is cut 
  ! ...      then reset variables
  ! ...      to their values at the end of the previous time step
!!$80 IF (tmcut.or.nswitch.gt.0) THEN
80 IF(tmcut) THEN
     ierr = .false.
     !****put this block at the end of the time step just after need for cut
     !*****    was determined
     ntimecut = ntimecut + 1
     ! ...      update the maximum number of N-R iterations for summary table
     ilnrimax = MAX(lnri, ilnrimax)
     ! ... Stop simulation if time step has been cut more than LTMCUTMX times
     IF (ntimecut > ltmcutmx) THEN
        WRITE (fustdout, 9010) ntimecut
        WRITE (fuclog, 9010) ntimecut
9010    FORMAT ('*** TIME STEP CUT', i3, ' TIMES - Simulation Stopped ***')
        ierr(110) = .TRUE.
        go to 900
     END IF
     ! ... *** if time step has been cut to a value < MINDELT - stop program ***
     ! ...          except for the first time step
     IF (delt < mindelt .AND. ltime > 1) THEN
        IF (iyr == 2) THEN
           WRITE(fuclog,9015) 'TIME STEP CUT TO  ',year(delt),' (yr)',  &
                'MINIMUM TIME STEP ',year(mindelt),' (yr)'
9015       FORMAT (tr5,a,1pg12.6,a/tr5,a,1pg12.6,a)
        ELSE
           WRITE(fuclog,9015) 'TIME STEP CUT TO  ',delt,' (s)',  &
                'MINIMUM TIME STEP ',mindelt,' (s)'
        END IF
        WRITE(fustdout,'(A)') '*** SIMULATION STOPPED - Time step too small ***'
        WRITE(fuclog,'(/A)') '*** SIMULATION STOPPED - Time step too small ***'
        ierr(111) = .TRUE.
        go to 900
     END IF
     IF (iyr == 2) THEN
        WRITE (fustdout,9020) ltime,'TIME STEP CUT,  new step length: ',year(delt),' (yr)'
        IF (ioptpr(3) >= 1) THEN
           WRITE (fuclog,9020) ltime,'TIME STEP CUT,  new step length: ',year(delt),' (yr)'
9020       FORMAT (/i4,tr4,a,1pg12.6,a)
        END IF
     ELSE
        WRITE (fustdout,9020) ltime,'TIME STEP CUT,  new step length: ',delt,' (s)'
        IF (ioptpr(3) >= 1) THEN
           WRITE (fuclog,9020) ltime,'TIME STEP CUT,  new step length: ',delt,' (s)'
        END IF
     END IF
     time = timeold + delt
     icall = 0
     iplot = 0
!!$     IF (prfreq < 0._kdp .AND. DBLE(ltime) >= DABS(prfreq)) iplot = 1
     ! ... If close to time for printout, move to time for printout
     ! ... If overshot time for printout, back up
     IF(ABS(time - prtimenxt) <= 0.1*delt .OR. time > prtimenxt) THEN
        deltoldt2 = delt
        delt = prtimenxt - timeold
        iplot = 1
     END IF
     ! ... *****This may end the time cut block to be moved
     tmcut = .FALSE.
     lnri = 0
     CALL reset
  END IF
  ! ... ******The logical end of the time cut block
  ! -------------------- NEWTON-RAPHSON ITERATION LOOP -------------------
  lnrichk1 = 4
  lnrichk2 = 6
  DO lnri=1,lnrimax
     lnritotl = lnritotl + 1
     IF (ioptpr(3) >= 1) THEN
        WRITE(fustdout,3001) 'N-R Iteration', lnri
3001    FORMAT(tr10,a,i3)
        WRITE(fuclog,3002) 'Newton-Raphson Iteration', lnri
3002    FORMAT(/tr10,a,i3)
     END IF
     DO  lll=1,npicmax
        ic = npic(lll)
        ! ...       Indoldnr(ic)=Ind(ic)
        poldnr(ic) = p(ic)
        holdnr(ic) = h(ic)
     END DO
     ! ... Assemble and solve linear NR matrix equation
     IF(slmeth == 1) THEN
        CALL ssor
     ELSEIF(slmeth == 2) THEN
        CALL gmres_asem_solve
        ! ... Check for error flags from Newton step with gmres
        DO ie=1,200
           IF(ierr(ie)) flag = 2
        END DO
        IF(flag > 0) RETURN
        IF (tmcut) THEN
           ! ... Cut time step if max number of GMRES iterations allowed is reached
           ! ...       without convergence
           delt = delt*0.5_kdp
           atmcut = 'G'
           GO TO 80
        END IF
     END IF
     ! ... Save maximum mass and energy residuals from previous N-R iterations
     IF (lnri >= 3) THEN
        resmold2 = resmold
        reseold2 = reseold
     END IF
     IF (lnri >= 2) THEN
        resmold = resmaxm
        reseold = resmaxe
     END IF
     ! ... Save the phase index from current N-R iteration
     DO  lll = 1, npicmax
        ic = npic(lll)
        indoldnr(ic) = ind(ic)
     END DO
     ! ... Use Cooley's algorithm for damping the change in dependent
     ! ...      variables; unconfined flow only
     IF(unconfined) THEN
        CALL damper
        !***..**restore arguments later F90
        !..          if(unconfined) call damper(wrelaxp,wrelaxh)
        IF(Ioptpr(3) >= 1) WRITE(FUCLOG, 9135) Ltime, Lnri, wrelaxp, wrelaxh
9135    FORMAT(I4,I3,' P under relaxation factor:',F6.2,'; H under relaxation factor:',F6.2) 
     ENDIF
     ! ... Update changes in pressure and enthalpy from N-R iteration
     DO  lll = 1, npicmax
        ic = npic(lll)
        p(ic) = poldnr(ic) + wrelaxp*dpn(ic)
        h(ic) = holdnr(ic) + wrelaxh*dhn(ic)
     END DO
     ! ...   check new P and H to make sure they are still within limits
     ! ...       if not, cut the time step
     DO  lll = 1, npicmax
        ic = npic(lll)
        !.. special patch for low pressures
        !..            IF (P(ic).GT.1.0D10 .OR. H(ic).GT.52.0D9 .OR.
        !..     &          P(ic).LT.0.5D6 .OR. H(ic).LT.0.01D9) THEN
        ! ... Test for (h,p) outside bounding box
        pmin = pp(1)
        pmax = pp(3)
        hmin = hh(1)
        hmax = hh(3)
        if(unconfined) pmin = 0._kdp
        IF (p(ic) < pmin .or. p(ic) > pmax .or. h(ic) < hmin .or. h(ic) > hmax) THEN
           ! ... point is outside of the bounding box of the table
           IF (ioptpr(3) >= 1) THEN
              CALL ictoijk(ic,i,j,k,nx,nz)
              WRITE (fustdout, 9025) ' *** P,H BEYOND LIMITS at IC =', ic,  &
                   '     I,J,K = ', i, j, k, ' ***', ' P= ', p(ic), '  H= ', h(ic)
9025          FORMAT (a, i6, a, 3i4, a/tr10,a , 1pg12.5, a, 1pg12.5)
              WRITE (fuclog, 9025) ' *** P,H BEYOND LIMITS at IC =', ic,  &
                   '     I,J,K = ', i, j, k, ' ***', ' P= ', p(ic), '  H= ', h(ic)
           END IF
!****                 ierr(71) = .TRUE.
           ! ...       increment timecut counter if it is the first time cut at this time step
           IF (.NOT.tmcut)  itmcntl(18) = itmcntl(18) + 1
           tmcut = .TRUE.          ! ... set flag for time step reduction
           atmcut = 'T'
        ELSE
           pmin = pp(2)
           hmax = hh(2)
           IF (.not.(p(ic) < pmin .or. p(ic) > pmax .or. h(ic) < hmin .or. h(ic) > hmax)) THEN
              ! ... point is outside of the table and inside the super-critical
              ! ...     compressed water phase
              IF (ioptpr(3) >= 1) THEN
                 CALL ictoijk(ic,i,j,k,nx,nz)
                 WRITE (fustdout, 9025) ' *** P,H BEYOND LIMITS at IC =', ic,  &
                      '     I,J,K = ', i, j, k, ' ***', ' P= ', p(ic), '  H= ', h(ic)
                 WRITE (fuclog, 9025) ' *** P,H BEYOND LIMITS at IC =', ic,  &
                      '     I,J,K = ', i, j, k, ' ***', ' P= ', p(ic), '  H= ', h(ic)
              END IF
!****                 ierr(71) = .TRUE.
              ! ...       increment timecut counter if it is the first time cut at this time step
              IF (.NOT.tmcut)  itmcntl(18) = itmcntl(18) + 1
              tmcut = .TRUE.          ! ... set flag for time step reduction
              atmcut = 'T'
           END IF
        END IF
        ! ...         ??? switch for checking for pressures greater than lithostatic, why???
        ! ... will never be done
        ioptlitho = .FALSE.
        ! ...   IOPTLITHO=.TRUE.
        ! ...         check for Pressure greater than lithostatic
        IF (ioptlitho .AND. p(ic) > plith(ic)) THEN
           CALL ictoijk(ic,i,j,k,nx,nz)      
           WRITE(fustdout,'(A,I4,A,3I4,A/1P,10X,A,G12.5,A,G12.5)')  &
                ' *** P > P lithostatic at IC =', ic,  &
                '     I,J,K = ', i, j, k, ' ***', ' P=  ', p(ic),  &
                '  PLITH=  ', plith(ic)
           WRITE(fuclog, '(A,I4,A,3I4,A/1P,10X,A,G12.5,A,G12.5)')  &
                ' *** P > P lithostatic at IC =', ic,  &
                '     I,J,K = ', i, j, k, ' ***', ' P=  ', p(ic),  &
                '  PLITH=  ', plith(ic)

           itmcntl(19) = itmcntl(19) + 1
           tmcut = .TRUE.
           atmcut = 'L'
        END IF
     END DO
     IF (tmcut) THEN
        delt = delt*0.1_kdp        ! ... Reduce time step to 10% of previous value
        GO TO 80
     END IF
     ! ... Determine the phase region (water, 2 phase, steam, super-critical fluid, or air-water)
     ! ...      for each cell
     CALL phreg(p,h)
     ! ... Determine if any phase changes occurred during current N-R iteration
     ! ...       or anytime since the previous time step
     CALL phasechg
     ! ... Check for error flags from phase changes
     DO ie=1,200
        IF(ierr(ie)) flag = 2
     END DO
     IF(flag > 0) RETURN
     IF(iphschgnr > 0) THEN
        lnrichk2 = MAX(lnri+3,6)
        lnrichk1 = lnrichk2-2
     ENDIF
     !???  if a phase change exceeded the saturation limit then cut the time step???
     ! ... If the time step has been cut, go back and repeat....this seems not possible here
     IF (tmcut) GO TO 80
     CALL tempdegc(p,h,tc)         ! ... Calculate temperature field for new P and H
     IF (lrxftreq) THEN
        CALL rockparm      ! ... Calculate permeability and transmissivity as f(Temperature)
        CALL tcalc     ! ... Calculate transmissivity and thermal conductivity for the cell faces
     END IF
     IF (lrxftimreq) THEN
        CALL permftime     ! ... Calculate permeability and transmissivity as f(Time)
        CALL tcalc     ! ... Calculate transmissivity and thermal conductivity for the cell faces
     END IF
     CALL enthrock          ! ... calculate enthalpy of the porous matrix
     CALL prpty(1)          ! ... calculate fluid property coefficients and derivatives
     CALL wellallo          ! ... determine the heat and mass allocations of source/sink terms
     IF (ioptupst) CALL upstre     ! ... determine the upstream node for each cell face
     CALL storativ          ! ... calculate mass and energy storativities
     ! ... calculate maximum changes in P and H and maximum percent changes
     CALL max_ddv_NR(pconvrgmx, pconvrgpc, hconvrgmx, hconvrgpc,  &
          icpdiffmax, ichdiffmax, icpdiffpc, ichdiffpc)
     CALL ictoijk(icpdiffmax,ipcon,jpcon,kpcon,nx,nz)
     CALL ictoijk(ichdiffmax,ihcon,jhcon,khcon,nx,nz)
     IF(ioptpr(3) >= 1) WRITE(fuclog,9035) ltime, lnri, pconvrgmx,  &
          ipcon, jpcon, kpcon, hconvrgmx, ihcon, jhcon, khcon
9035 FORMAT (i4,i3,' Max Diff: ',1pg11.3,' in P at ',3I4,'; ',  &
          1pg11.3,' in H at ',3I4)
     CALL ictoijk(icpdiffpc,ippc,jppc,kppc,nx,nz)
     CALL ictoijk(ichdiffpc,ihpc,jhpc,khpc,nx,nz)
     !  
     IF(ioptpr(3) >= 1) THEN
        IF(unconfined) THEN
           WRITE(fuclog,9140) ltime, lnri, pconvrgmx,  &
                ipcon, jpcon, kpcon, hconvrgpc, ihpc, jhpc, khpc
9140       FORMAT (i4,i3,' Max Diff: ',1pg11.3,' in P at ',3I4,';',  &
                1pg11.3,'% in H at ',3I4)
        ELSE
           WRITE(fuclog,9040) ltime, lnri, pconvrgpc,  &
                ippc, jppc, kppc, hconvrgpc, ihpc, jhpc, khpc
9040       FORMAT (i4,i3,' Max %diff: ',1pg11.3,' in P at ',3I4,';',  &
                1pg11.3,' in H at ',3I4)
        END IF
     END IF
     CALL resid          ! ... calculate mass and energy residuals
     IF (ioptpr(3) >= 1) WRITE (fuclog,9030) ltime, lnri, resmaxm,  &
          imres, jmres, kmres, resmaxe, ieres, jeres, keres
     ! ... ------------- Test for N-R iteration convergence -------------
     ! ...   Determine if the residuals after N-R iteration meet convergence criteria
     ! ...       Resmas - input value for maximum acceptable residual mass
     ! ...       Reseng - input value for maximum acceptable residual energy
     ! ... Convergence criteria: 3 possibilities
     ! ...    1)  if residuals are sufficiently small
     ! ...    2)  on/after 5 NRI, if absolute changes in P&H are sufficiently small
     ! ...    3)  on/after 4 NRI, if % change in P&H is sufficiently small
     ! ...         and residuals are within a factor of 100 of their convergence criteria
     ! ...         and a phase change has occurred
     !
     ! ...       determine which of the above criteria were met for convergence
     IF (ABS(resmaxm) <= resmas .AND. ABS(resmaxe) <= reseng) THEN
        lconvrg = 1
     ELSE IF (lnri >= 5 .AND.  &
          ABS(pconvrgmx) <= 1._kdp .AND. ABS(hconvrgmx) <= 1.e2_kdp) THEN
        lconvrg = 2
     ELSE IF(.NOT.unconfined .AND. lnri >= 4 .AND. iphschglt >= 1 .AND.  &
          ABS(pconvrgpc) <= pchgmxnr .AND. ABS(hconvrgpc) <= hchgmxnr .AND.  &
          ABS(resmaxm) <= (100._kdp*resmas) .AND.  &
          ABS(resmaxe) <= (100._kdp*reseng)) THEN
        lconvrg = 3
     ELSE IF(unconfined .AND. lnri >= 4 .AND. iphschglt >= 1 .AND.  &
          ABS(pconvrgmx) <= pchgmxnr .AND. ABS(hconvrgpc) <= hchgmxnr .AND.  &
          ABS(resmaxm) <= (100._kdp*resmas) .AND.  &
          ABS(resmaxe) <= (100._kdp*reseng)) THEN
        lconvrg = 3
     ELSE
        ! ... No convergence
        lconvrg = 0
        ! **** flag = ??
     END IF
     IF(lconvrg > 0) THEN
        ! ...       convergence achieved.
        ! ...       determine the max amount of water saturation change and
        ! ...      the max % change in P and H over the current time step
        CALL maxchg(wschgocr, pchgocr, hchgocr, iwschg, jwschg, kwschg,  &
             ipchg, jpchg, kpchg, ihchg, jhchg, khchg, atmcut)
        ! ... MAXCHG adjusts the time step based on maximum dependent variable changes
        IF (tmcut) GO TO 80
        ! ... Keep a record of which convergence criteria were met
        CALL converge_criterion(lconvrg, aconvrg, uconvrg)
        ! ... If no time step cut and no seepage switches, then
        ! ...      break out of N-R iteration loop.
        GO TO 140        ! ... This is the successful exit point
     END IF
     ! ... If we are here, then residual convergence criteria are not met.
     ! ... Determine if a smaller time step should be used if it is the 9th 
     ! ...     N-R iteration. Check the max change in water sat'n
     ! ...      and the max % change in P and H since the last time step
     ! ...      Or the max abs change in P if unconfined flow
     IF (lnri > 8) THEN
        CALL maxchg(dum1, dum2, dum3, idum4, idum5, idum6,  &
             idum7, idum8, idum9, iduma, idumb, idumc, atmcut)
        ! ... MAXCHG adjusts the time step based on maximum dependent variable changes
        IF (tmcut) GO TO 80
     END IF
     ! ...   cut time step if more than IPHSCHGMX nodes changed phase
     ! ...        during the current N-R iteration
     IF (iphschgnr > iphschgmx) THEN
        IF (ioptpr(3) >= 1) THEN
           WRITE (fustdout,'(tr8,i3,2a)') iphschgnr,  &
                ' Nodes Changed Phase during one N-R iteration; ',  &
                'EXCEEDS MAXIMUM ALLOWED'
           WRITE (fuclog,'(tr8,i3,2a)') iphschgnr,  &
                ' Nodes Changed Phase during one N-R iteration; ',  &
                'EXCEEDS MAXIMUM ALLOWED'
        END IF
        tmcut = .TRUE.
        itmcntl(14) = itmcntl(14) + 1
        delt = delt*0.5_kdp        ! ... reduce time step by 0.5 factor
        atmcut = 'C'
        GO TO 80
     END IF
     ! ...   cut time step if more than IPHSCHGMX+2? nodes changed phase
     ! ...        during the current time step
     IF (iphschglt > iphschgmx+1) THEN
        IF (ioptpr(3) >= 1) THEN
           WRITE (fustdout,'(tr8,i3,a)') iphschglt,  &
                ' Nodes Changed Phase over one time step; EXCEEDS MAXIMUM ALLOWED'
           WRITE (fuclog,'(tr8,i3,a)') iphschglt,  &
                ' Nodes Changed Phase over one time step; EXCEEDS MAXIMUM ALLOWED'
        END IF
        tmcut = .TRUE.
        itmcntl(15) = itmcntl(15) + 1
        delt = delt*0.5_kdp        ! ... reduce time step by 0.5 factor
        atmcut = 'C'
        GO TO 80
     END IF
     ! ...   cut time step if max number of N-R iterations allowed is reached
     ! ...        without convergence
     IF (lnri >= lnrimax) THEN
        IF (ioptpr(3) >= 1) THEN
           WRITE (fustdout,'(2A)') '*** MAXIMUM NUMBER of NEWTON-RAPHSON ',  &
                'ITERATIONS REACHED WITHOUT CONVERGENCE ***'
           WRITE (fuclog,'(2A)') '*** MAXIMUM NUMBER of NEWTON-RAPHSON ',  &
                'ITERATIONS REACHED WITHOUT CONVERGENCE ***'
        END IF
        tmcut = .TRUE.
        itmcntl(16) = itmcntl(16) + 1
        delt = delt*0.5_kdp
        atmcut = 'I'
        GO TO 80
     END IF
     !******remove this test for monotonic decreases ****
!!$     ! ...    For NR iterations #4 and #5, check if mass or energy balance errors
!!$     ! ...      have increased since the two previous  N-R iterations.
!!$     ! ...    Beginning with NR iteration #6, check if residuals have increased
!!$     ! ...      relative to the residuals two iterations prior.
!!$     ! ... Residuals allowed to increase if there was a phase change
!!$     IF (Lnri == lnrichk1 .or. Lnri == lnrichk1+1) THEN
!!$        IF(ABS(Resmaxm) > ABS(resmold) .AND. ABS(Resmaxe) > ABS(reseold) .AND. &
!!$             ABS(Resmaxm) > ABS(resmold2) .AND. ABS(Resmaxe) > ABS(reseold2)) then
!!$           IF (ioptpr(3) >= 1) THEN
!!$              WRITE (fustdout, '(A)')  &
!!$                   '*** NEWTON-RAPHSON ITERATION NOT CONVERGING ***'
!!$              WRITE (fuclog, '(A)') '*** NEWTON-RAPHSON ITERATION NOT CONVERGING ***'
!!$           END IF
!!$           tmcut = .true.
!!$           itmcntl(17) = itmcntl(17) + 1
!!$           delt = delt*0.5_kdp
!!$           atmcut = 'N'
!!$           GO TO 80
!!$        END IF
!!$     ELSEIF(Lnri >= lnrichk2 .and. ABS(Resmaxm) > ABS(resmold2) .AND. &
!!$          ABS(Resmaxe).GT.ABS(reseold2)) THEN
!!$        IF (ioptpr(3) >= 1) THEN
!!$           WRITE (fustdout, '(A)')  &
!!$                '*** NEWTON-RAPHSON ITERATION NOT CONVERGING ***'
!!$           WRITE (fuclog, '(A)') '*** NEWTON-RAPHSON ITERATION NOT CONVERGING ***'
!!$        END IF
!!$        tmcut = .TRUE.
!!$        itmcntl(17) = itmcntl(17) + 1
!!$        delt = delt*0.5_kdp
!!$        atmcut = 'N'
!!$        GOTO 80
!!$     ENDIF
  END DO
  ! ... -------------- END OF NEWTON-RAPHSON ITERATION LOOP --------------
140 CONTINUE
  ! ... Save maximum residual values for summary table
  IF (ABS(resmaxm) > ABS(resmaxmall)) THEN
     resmaxmall = resmaxm
     ltimresmaxm = ltime
  END IF
  IF (ABS(resmaxe) > ABS(resmaxeall)) THEN
     resmaxeall = resmaxe
     ltimresmaxe = ltime
  END IF
  ! ...       write convergence statement
  !*****needs cleaning up
  IF (ioptpr(3) == 0) THEN
     IF (iyr == 2) THEN
        IF (ntimecut > 0) THEN
           WRITE (fuclog, 9045) ltime, ' Time Step Length (yr): ', year(delt),  &
                year(time), lnri, atmcntl, aconvrg, iphschglt, ';', ntimecut, atmcut
        ELSE IF (iphschglt > 0) THEN
           WRITE (fuclog, 9045) ltime, ' Time Step Length (yr): ', year(delt),  &
                year(time), lnri, atmcntl, aconvrg, iphschglt, ';'
        ELSE
           WRITE (fuclog, 9045) ltime, ' Time Step Length (yr): ', year(delt),  &
                year(time), lnri, atmcntl, aconvrg
        END IF
     ELSE IF(iyr == 1) THEN
        IF (ntimecut > 0) THEN
           WRITE (fuclog, 9045) ltime, ' Time Step Length (s): ', delt, time,  &
                lnri, atmcntl, aconvrg, iphschglt, ';', ntimecut, atmcut
        ELSE IF (iphschglt > 0) THEN
           WRITE (fuclog, 9045) ltime, ' Time Step Length (s): ', delt, time,  &
                lnri, atmcntl, aconvrg, iphschglt, ';'
        ELSE
           WRITE (fuclog, 9045) ltime, ' Time Step Length (s): ', delt, time,  &
                lnri, atmcntl, aconvrg
        END IF
     ENDIF
9045 FORMAT (i5,a,1pg9.3,' Time: ',1pg10.4,i3,' NR iter ', a1,'-cntl ',  &
          a1, '-cnvrg;', :,i2,' PhsChg',a1,i2,' cut-',a1)
  END IF
  IF (ioptpr(3) >= 1) THEN
     WRITE (fustdout,'(I4,4X,A,A,A,I3,A,I3,A)') ltime,  &
          'Convergence: criterion met- ',uconvrg,';',lnri,  &
          ' NR iters;',iphschglt,' nodes changed to/from 2-phase'
     WRITE (fuclog, '(I4,4X,A,A,A,I3,A,I3,A)') ltime,  &
          'Convergence: criterion met- ',uconvrg,';',lnri,  &
          ' NR iters;',iphschglt,' nodes changed to/from 2-phase'
     IF (ioptpr(3) == 1) THEN
        IF(unconfined) THEN
           WRITE (fuclog,9150) ltime, pchgocr, hchgocr
9150       FORMAT (i4, 4X, 'Max change in P = ',1pg11.4,  &
                ';  Max % change in H = ',0pf10.3)
        ELSE
           WRITE (fuclog,9050) ltime, pchgocr, hchgocr
9050       FORMAT (i4, 4X, 'Max % change in P = ',0pf10.3,  &
                ';  Max % change in H = ',0pf10.3)
        END IF
     ELSEIF(ioptpr(3) == 2) THEN
        icp = (kpchg-1)*nx + ipchg + (jpchg-1)*nx*nz
        ich = (khchg-1)*nx + ihchg + (jhchg-1)*nx*nz
        IF(unconfined) THEN
           WRITE(fustdout,9155) ltime, p(icp) - poldt(icp), ipchg,  &
                jpchg, kpchg, ltime, hchgocr, h(ich) - holdt(ich), ihchg, jhchg, khchg
           WRITE(fuclog,9155) ltime, p(icp) - poldt(icp), ipchg,  &
                jpchg, kpchg, ltime, hchgocr, h(ich) - holdt(ich), ihchg, jhchg, khchg
9155       FORMAT(i4,tr4,'Max delta P       = ',1pg11.4,' at I,J,K ',3I4/i4,tr4,  &
                'Max % change in H = ',0pf10.3,'; delta H = ',1pg11.4,' at I,J,K ',3I4)
        ELSE
           WRITE(fustdout,9055) ltime, pchgocr, p(icp) - poldt(icp), ipchg,  &
                jpchg, kpchg, ltime, hchgocr, h(ich) - holdt(ich), ihchg, jhchg, khchg
           WRITE(fuclog,9055) ltime, pchgocr, p(icp) - poldt(icp), ipchg,  &
                jpchg, kpchg, ltime, hchgocr, h(ich) - holdt(ich), ihchg, jhchg, khchg
9055       FORMAT(i4,tr4,'Max % change in P = ',f10.3,'; delta P = ',  &
                1pg11.4,' at I,J,K ',3I4/i4,tr4,  &
                'Max % change in H = ',0pf10.3,'; delta H = ',1pg11.4,' at I,J,K ',3I4)
        END IF
     END IF
  END IF
  ! ... calculate amount of vapor, liquid, and super-critical fluid in system
  CALL steam_mass(vpir, lqir, scrfir)
  IF (wschgocr > 0._kdp .OR. vpir > 0._kdp) THEN
     IF (ioptpr(3) == 1) THEN
        WRITE(fuclog, 9060) ltime, wschgocr, vpir
     ELSEIF (ioptpr(3) == 2) THEN
        WRITE(fuclog, 9065) ltime, wschgocr, iwschg, jwschg, kwschg, ltime, vpir
     END IF
  END IF
  ! ...       calculate the mass and heat balance for this time step
  CALL balance
  ! ... Print current time step length and time on screen
  IF (iyr == 2) THEN
     WRITE (fustdout,9070) ltime,'Time Step Length: ',year(delt),' (yr)'
9070 FORMAT (i4,tr25,a,1pg12.6,a)
     WRITE (fustdout,9075) 'Step Number',ltime,  &
          ' completed; Simulation time .................... ',  &
          year(time),' (yr)'
9075 FORMAT (a,i5,a,1pg12.6,a)
     IF (ioptpr(3) >= 1) THEN
        WRITE (fuclog,9075) 'Step Number',ltime,  &
             ' completed; Simulation time .................... ',  &
             year(time),' (yr)'
     END IF
  ELSE
     WRITE (fustdout,9070) ltime,'Time Step Length: ',delt,' (s)'
     WRITE (fustdout,9075) 'Step Number',ltime,  &
          ' completed; Simulation time .................... ',  &
          time,' (s)'
     IF (ioptpr(3) >= 1) THEN
        WRITE (fuclog,9075) 'Step Number',ltime,  &
             ' completed; Simulation time .................... ',  &
             time,' (s)'
     END IF
  END IF
  ! ... Save information for simulation summary
  ! ...    max and min time step lengths
  deltmax = MAX(delt,deltmax)
  deltmin = MIN(delt,deltmin)
  ! ...    maximum number of N-R iterations
  ilnrimax = MAX(lnri,ilnrimax)
  ! ... total number of phase changes during simulation
  iphschgtot = iphschgtot + iphschglt
  ! ... Print data at selected time planes or at time of change in the
  ! ...      input parameters or if error abort
  ! ... Calculate the water and steam velocity for printing if necessary
  IF (ABS(print_vel) > 0._kdp .OR. ABS(print_plotvector) > 0._kdp .OR.  &
       ABS(print_plotscalar) > 0._kdp .OR.  &
       ioptprta(9) >= 1 .OR. ioptprta(10) >= 1 .OR.  &
       ioptprta(11) >= 1 .OR. ioptprta(15) >= 1 .OR.  &
       ioptprta(16) >= 1 .OR. ioptprta(17) >= 1 .OR. ioptprta(27) >=1) THEN
     CALL velocity
  END IF
  ! ... Write tabular information to output files
  CALL pdata(flag)
  IF(ABS(print_plotscalar) > 0._kdp .OR. ABS(print_plotvector) > 0._kdp) THEN
     SELECT CASE (ioptpr(1))
     CASE (1)        ! ... write Explorer plotfile
        CALL plotexpl
     CASE (2,3)      ! ... write gnuplot plotfile
        CALL plotgnu
     CASE (4,5)      ! ... write HTView IDL plotfile
        CALL plotidl
     CASE (6)        ! ... Write columnar xyz plotfile
        CALL plotxyz
     CASE (8,9)      ! ... write B2plotfile
        CALL plotb2
     END SELECT
  END IF
  IF(ABS(print_dump) > 0._kdp) THEN
     ! ... Print restart information to restart file
     CALL print_control(print_dump,time,ltime,tchg,timprdump,lprint)
     IF(flag > 0) lprint = .TRUE.          ! ... force print if an error condition
     IF(lprint) CALL pickup
  END IF
  ! ... Set the time for next printout
  prtimenxt = MIN(tchg, timprp, timprh, timprt, timprsat, timprd, timprv, timprpot,  &
       timprvel, timprbcf, timprsrc, timprprop, timprpor, timprperm, timprbal, timprdump,  &
       timprresid, timprdimno, timprpscal, timprpvec, timprtemp)
  nswitch = 0
  IF(nseep > 0) THEN
     ! ... Switch seepage b.c. nodes from active to passive as necessary
     DO  j = 1,ny
        DO  k = 1,nz
           DO  i = 1,nx
              ij = (j-1)*nx + i
              ik = (k-1)*nx + i
              ! ... Select seepage nodes or combination seepage/precipitation nodes
              IF (ibc(ik,j) == 30 .OR. ibc(ik,j) == 50) THEN
                 ic = (k-1)*nx + i + (j-1)*nxz
                 IF(seeping(ik,j)) THEN
                    ! ... Switch active seepage nodes based on net flux
                    CALL gov_seep(gveqnm, gveqne, i, j, k)
                    ! ... Add in the capacitance term
                    ! ...     No source included for seeping cell
                    gveqnm = gveqnm - (xm(ic)-xmoldt(ic))/delt
!!$                       qseep(ic) = -gveqnm     ! ... Net flow in is seepage outflow
!!$                          gveqne = gveqne - (en(ic)-enoldt(ic))/delt
                    ! ... Convert to flux and get correct sign (out is negative)
                    qfluxsp = -gveqnm/land_surf_area(ij)
!!$                          qhfluxsp = -gveqne/land_surf_area(ij)
                    ! ... Net flux into cell means switch to noflux b.c. or recharge flux
                    ! ...     b.c.; passive seepage node
                    IF (qfluxsp > 0._kdp) THEN
                       seeping(ik,j) = .FALSE.
                       nswitch = nswitch + 1
                    ENDIF
                 ELSEIF(.NOT.seeping(ik,j)) THEN
                    ! ...      Switch passive seepage nodes based on fluid pressure
                    IF(p(ic) >= patm) THEN
                       p(ic) = patm
                       seeping(ik,j) = .TRUE.
                       nswitch = nswitch + 1
                    ENDIF
                 END IF
              ENDIF
           END DO
        END DO
     END DO
  ENDIF
!!$      *** following is suspended for the present
!!$        ! ... If any seepage cells switched type, repeat the time step
!!$        IF (nswitch > 0) then
!!$           IF(Ioptpr(3) >= 1) WRITE(FUCLOG, 9235) Ltime, Lnri, &
!!$                'Repeat time step from seepage switch'
!!$9235       FORMAT (I4,I4,tr2,a) 
!!$           GOTO 80
!!$        ENDIF
  ! ... Save old time step (DELTOLDT)
  ! ... Set new time step based on the max change in P, H, and water saturation
  deltoldt2 = deltoldt
  deltoldt = delt
  awschgocr = MAX(ABS(wschgocr),1.e-25_kdp)
  apchgocr = MAX(ABS(pchgocr),1.e-25_kdp)
  ahchgocr = MAX(ABS(hchgocr),1.e-25_kdp)
  ! ... Factor of 0.5 has been removed. It made the targets half of the values
  ! ...      specified
  fac = MIN(tmincmx, pchgmxts/apchgocr, hchgmxts/ahchgocr, wschgmx/awschgocr)
  ! ... Increase time step only if it is less than 0.9 of desired value
  ! ... Accept any decrease in time step
!    due to machine truncation error, fac could be slight less than 1.
  IF(fac > 0.99999_kdp .AND. fac < 1.11_kdp) fac = 1._kdp
  delt = deltoldt*fac
  ! ... Record which criterion determined the new time step
  IF (fac == 1._kdp) THEN
     ltmcntl = 0
  ELSEIF (fac == wschgmx/awschgocr) THEN
     ltmcntl = 1
  ELSEIF (fac == pchgmxts/apchgocr) THEN
     ltmcntl = 2
  ELSEIF (fac == hchgmxts/ahchgocr) THEN
     ltmcntl = 3
  ELSEIF (fac == tmincmx) THEN
     ltmcntl = 4
  END IF
  ! ... If current time step was set because it was time to print, restore the previous
  ! ...      time step if it is greater than the
  ! ...      new time step determined by saturation, pressure, or enthalpy change
!!$ IF (iplot == 1 .AND. icall == 0 .AND. deltoldt2*tmincmx >= delt .AND. prfreq > 0._kdp) THEN
  IF(iplot == 1 .AND. icall == 0 .AND. deltoldt2 >= delt) THEN
     delt = deltoldt2
     ltmcntl = 5
  END IF
  ! ... If time for new period, read in new time step, source, and parameter data
  IF (icall == 1) THEN
     delt = deltoldt2
     CALL sink(nppn,flag)
     ! ... SINK may read the termination flag and prepare to end the simulation
     IF(flag > 0) RETURN
     ! ... DELT may or may not have been reset in SINK
     IF (delt == deltoldt2) THEN
        ltmcntl = 9
     ELSE
        ltmcntl = 6
     END IF
     maxdelt = tchg     ! ... *** temporary 
  END IF
  iplot = 0
  ! ... If close to time for printout, move to time for printout
  ! ... If overshot time for printout, back up
  IF(ABS(time+delt-prtimenxt) <= 0.1*delt .OR. time+delt > prtimenxt) THEN
     deltoldt2 = delt
     delt = prtimenxt - time
     iplot = 1
     ltmcntl = 8
  END IF
!!$  ! ...          if printing occurred this time step set next time to print
!!$  IF (iplot == 1 .AND. prfreq > 0._kdp) THEN
!!$     prtime = (1._kdp + DINT(time/prfreq))*prfreq
!!$     ! ... If time of next printout is within 10% of the print interval to
!!$     ! ...      the time for read in of 
!!$     ! ...      new data, move printout time to read in time.
!!$     IF(ABS(tchg-prtime) < 0.1_kdp*prfreq) prtime = tchg
!!$  END IF
!!$  ! ... If simulation time is about to exceed the time to print data,
!!$  ! ...      adjust the time step to hit the print time.
!!$  ! ...      Only for print interval output control, not time step
!!$  ! ...           output control.
!!$  iplot = 0
!!$  IF((time+delt) >= prtime) THEN
!!$     iplot = 1
!!$     ltmcntl = 8
!!$     delt = prtime - time
!!$  ENDIF
  icall = 0
  ! ... If close to time for change, move to time for change or
  ! ...     if overshot time for change, back up
  IF(ABS(time+delt-tchg) <= 0.2*delt .OR. time+delt > tchg) THEN
     deltoldt2 = delt
     delt = tchg - time
     icall = 1
     ltmcntl = 7
  END IF
!!$  ! ... If simulation time is about to exceed the time to read data,
!!$  ! ...      adjust the time step to hit the data read-in time.
!!$  icall = 0
!!$  IF(time+delt >= tchg) THEN
!!$     delt = tchg - time
!!$     icall = 1
!!$     ltmcntl = 7
!!$  END IF
  ! ... Enforce maximum time step length
  IF(delt > maxdelt) THEN
     delt = maxdelt
     ltmcntl = 10
  END IF
  ! ... ####### End of a time step simulation  #######
  ! ... Check for error flags
900 CONTINUE
  DO ie=1,200
     IF(ierr(ie)) flag = 2
  END DO
  !*****maybe should check for error flags after every NR iteration
  !
9030 FORMAT (i4, i3,' Max Resid:',1PE10.2,' g/s-cm^3 at',3I4, ';',  &
       1PE11.2, ' erg/s-cm^3 at',3I4)
9060 FORMAT (i4, 4X,'Max Saturation change: ',f7.4,  &
       ';  Steam in region: ',1pg10.4,' (g)')
9065 FORMAT (i4, 4X,'Max Saturation change: ',f8.5,' at I,J,K ',  &
       3I4/i4,tr4,'Steam in region: ',1pg11.5,' (g)')
END SUBROUTINE time_step_ht
