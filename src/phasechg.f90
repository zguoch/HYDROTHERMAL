SUBROUTINE phasechg
  ! ... PURPOSE:  To check if any nodes have changed phase during a
  ! ...      NR-iteration
  ! ...      if they have, adjust the pressure and
  ! ...      enthalpy to near the two-phase curve.  This adjustment
  ! ...      helps with NR convergence.
  ! ...      This is not done for water/air-water phase changes.
  ! ... METHOD:  If a node changes phase on the first NR iteration, adjust
  ! ...      the P & H to a point just across the phase boundary.  This is
  ! ...      accomplished by stepping along the line between the old P & H
  ! ...      and the new P & H until the point is just inside the boundary.
  ! ...      If a node changes phase on the second NR iteration, suppress
  ! ...      the phase change by adjusting the P & H to a point next to,
  ! ...      but not across, the phase boundary.
  USE machine_constants, ONLY: kdp
  USE bc, ONLY: ibc
  USE f_units
  USE control
  USE mesh
  USE variables
  IMPLICIT NONE
  CHARACTER(LEN=14) :: newphase, oldphase
  INTEGER :: i, ic, ik, j, k, ml
  REAL(KIND=kdp) :: ainc, aval, dum1, hs, htemp, hw, ptemp
  !
  ! ...   Iphschgnr - number of nodes changing phase this NR-iteration
  ! ...   Iphschglt - number of nodes changing phase this time step
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.7 $//$Date: 2007/01/16 21:24:56 $'
  !     ------------------------------------------------------------------
  ! ... Set switch to determine if any cells had phase changes over this N-R
  ! ...     iteration
  iphschglt = 0
  iphschgnr = 0
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF (ibc(ik,j) == -1) CYCLE        ! ... Skip inactive cells
! ***** unnecessary test, done in phreg
!!$     IF (ind(ic) < 1 .OR. ind(ic) > 5) THEN     ! ... skip node if ind(ic) is invalid index
!!$        WRITE(fustdout, '(A,3I4/10X,A,I1)')  &
!!$             '*** STOP ***   Phase Index not valid for node I,J,K ',  &
!!$             i,j,k,' IND=  ',ind(ic)
!!$        WRITE(fuclog, '(A,3I4/10X,A,I1)')  &
!!$             '*** STOP ***   Phase Index not valid for node I,J,K ',  &
!!$             i,j,k,' IND=  ',ind(ic)
!!$        ierr(102) = .TRUE.
!!$        CYCLE
!!$     END IF
     ! ... Count the number of phase changes for the current time step
     ! ... Count only changes involving 2-phase region or air-water cells
     IF (ind(ic) /= indoldt(ic) .AND.  &
          ((indoldt(ic) == 2 .OR. ind(ic) == 2) .OR.  &
          (indoldt(ic) == 5 .OR. ind(ic) == 5))) iphschglt = iphschglt + 1
     IF (ind(ic) /= indoldnr(ic)) THEN     
        ! ... Count and print message if cell has changed phase
        iphschgnr = iphschgnr + 1
        SELECT CASE (indoldnr(ic))
        CASE (1)
           oldphase = 'Water'
        CASE (2)
           oldphase = 'Two-Phase'
        CASE (3)
           oldphase = 'Steam'
        CASE (4)
           oldphase = 'Super Critical'
        CASE (5)
           oldphase = 'Air-Water'
        END SELECT
        SELECT CASE (ind(ic))
        CASE (1)
           newphase = 'Water'
        CASE (2)
           newphase = 'Two-Phase'
        CASE (3)
           newphase = 'Steam'
        CASE (4)
           newphase = 'Super Critical'
        CASE (5)
           newphase = 'Air-Water'
        END SELECT
        ! ...     if 2-phase region not involved in "phase change" then don't
        ! ...       count it as a phase change. i.e. ignore changes between the
        ! ...       super-critical region and the water region, and between the
        ! ...       super-critical region and the steam region.
        ! ...       or if air-water region not involved, do not count as phase change
        IF ((indoldnr(ic) /= 2 .AND. ind(ic) /= 2) .AND.  &
             (indoldnr(ic) /= 5 .AND. ind(ic) /= 5)) THEN
           iphschgnr = iphschgnr - 1
           IF (ioptpr(3) >= 1) THEN
              WRITE(fustdout,'(I4,I4,A,3I4)') ltime, lnri,  &
                   '     Node Changed from '//TRIM(oldphase)//' to '//  &
                   TRIM(newphase)//'     at I,J,K= ', i, j, k
              WRITE(fuclog,'(I4,I4,A,3I4)') ltime, lnri,  &
                   '     Node Changed from '//TRIM(oldphase)//' to '//  &
                   TRIM(newphase)//'     at I,J,K= ', i, j, k
              WRITE(fuclog,9010) '  P & H of node: ',p(ic),h(ic),' '
9010          FORMAT(tr10,a,2(1pe13.6),tr1,a)
           END IF
           CYCLE
        ELSE IF (ioptpr(3) >= 1) THEN
           ! ... Print phase chages involving the 2-phase region and air-water region
           WRITE(fustdout,'(I4,I4,A,3I4)') ltime, lnri,  &
                '     Node Changed from '//TRIM(oldphase)//' to '//  &
                TRIM(newphase)//'     at I,J,K= ', i, j, k
           WRITE(fuclog,'(I4,I4,A,3I4)') ltime, lnri,  &
                '     Node Changed from '//TRIM(oldphase)//' to '//  &
                TRIM(newphase)//'     at I,J,K= ', i, j, k
           WRITE(fuclog,9010) '  P & H of node: ',p(ic),h(ic),' '
        END IF
        IF(indoldt(ic) == 5 .OR. ind(ic) == 5) CYCLE     ! ... No adjustments for air-water
        IF (lnri >= 2) THEN
           ! ... Write information and go to next node since this is not the first NR iteration
           IF (ioptpr(3) >= 1) THEN
              IF (p(ic) < 220.55e6_kdp) THEN
                 CALL phsbndy1(p(ic),hw,hs)
                 IF (indoldnr(ic) == 1 .OR. ind(ic) == 1) THEN
                    WRITE(fuclog,9005) '  P & H of node: ',p(ic),  &
                         h(ic), '; At this P, HW=', hw
                 ELSE
                    WRITE(fuclog,9005) '  P & H of node: ',p(ic),  &
                         h(ic), '; At this P, HS=', hs
9005                FORMAT(tr10,a,2(1pe13.6),tr1,a,1pe13.6)
                 END IF
              END IF
           END IF
           CYCLE
        END IF
        ! ... ################### FIRST NR ITERATION #########################
        ! ...       write P & H of node at end of previous time step and at this NR iteration
        ! ...       also write the saturated H for water or steam at each P
        IF (ioptpr(3) >= 1) THEN
           IF (poldt(ic) < 220.55e6_kdp) THEN
              CALL phsbndy1(poldt(ic), hw, hs)
              IF (indoldnr(ic) == 1 .OR. ind(ic) == 1) THEN
                 WRITE(fuclog, 9005) 'P&H previous time step',poldt(ic),  &
                      holdt(ic),'; At this P, HW=', hw
              ELSE
                 WRITE(fuclog, 9005) 'P&H previous time step',poldt(ic),  &
                      holdt(ic),'; At this P, HS=', hs
              END IF
           END IF
           IF (p(ic) < 220.55e6_kdp) THEN
              CALL phsbndy1(p(ic), hw, hs)
              IF (indoldnr(ic) == 1 .OR. ind(ic) == 1) THEN
                 WRITE(fuclog, 9005) 'P&H will be adjusted from ',p(ic),h(ic),  &
                      '; At this P, HW=', hw
              ELSE
                 WRITE(fuclog, 9005) 'P&H will be adjusted from ',p(ic),h(ic),  &
                      '; At this P, HS=', hs
              END IF
           END IF
        END IF
        ! ...      set values for the beginning point and increment amount
        ! ...        for the P & H adjustment
        aval = 0.50_kdp
        ainc = 0.01_kdp
        IF (indoldnr(ic) == 1 .AND. ind(ic) == 2) THEN
           ! ... ****************  PHASE CHANGE:  WATER to Two-PHASE ****************
           ! ...      Adjust the pressure and enthalpy of any node that has gone from
           ! ...        compressed water at the last N-R iteration to two phase region
           ! ...        at this iteration
           ! ...        (and was in the compressed water region last time step)
           ! ...        to a P and H just inside the two phase/water boundary.
           ! ...        Step along the line between the old P & H and the new
           ! ...        P & H until the point is just inside the boundary.
           ptemp = poldt(ic)
           htemp = holdt(ic)
           IF (ptemp >= 220.55e6_kdp) THEN
              hw = 20.86e9_kdp
           ELSE
              CALL phsbndy1(ptemp,hw,dum1)
           END IF
           aval = 0.0_kdp
           DO WHILE (htemp < hw .OR. ptemp >= 220.55e6_kdp)
              aval = aval + ainc
              ptemp = aval*p(ic) + (1.0_kdp-aval)*poldt(ic)
              htemp = aval*h(ic) + (1.0_kdp-aval)*holdt(ic)
              IF (ptemp >= 220.55e6_kdp) THEN
                 hw = 20.86e9_kdp
              ELSE
                 CALL phsbndy1(ptemp,hw,dum1)
              END IF
           END DO
           p(ic) = ptemp
           h(ic) = htemp
        ELSEIF (indoldnr(ic) == 2 .AND. ind(ic) == 1) THEN
           ! ... ****************  PHASE CHANGE:  Two-PHASE to WATER ****************
           ! ...      Adjust the pressure and enthalpy of any node that has gone from
           ! ...        two phase region at the last NR iteration to compressed water at 
           ! ...        this iteration
           ! ...        to a P and H just outside the two phase/water boundary.
           ! ...        Step along the line between the old P & H and the new
           ! ...        P & H until the point is just outside the boundary.
           ptemp = poldt(ic)
           htemp = holdt(ic)
           IF (ptemp >= 220.55e6_kdp) THEN
              hw = 20.86e9_kdp
           ELSE
              CALL phsbndy1(ptemp, hw, dum1)
           END IF
           aval = 0.0_kdp
           DO                              
              IF(htemp < hw) EXIT
              aval = aval + ainc
              ptemp = aval*p(ic) + (1.0_kdp-aval)*poldt(ic)
              htemp = aval*h(ic) + (1.0_kdp-aval)*holdt(ic)
              IF (ptemp >= 220.55e6_kdp) THEN
                 hw = 20.86e9_kdp
              ELSE
                 CALL phsbndy1(ptemp, hw, dum1)
              END IF
           END DO
           p(ic) = ptemp
           h(ic) = htemp
        ELSEIF (indoldnr(ic) == 2 .AND. ind(ic) >= 3) THEN
           ! ... ****  PHASE CHANGE:  Two-PHASE to SUPERHEATED STEAM or SUPER CRITICAL ****
           ! ...      Adjust the pressure and enthalpy of any node that has gone from
           ! ...        the two phase region at the last N-R iteration to either superheated
           ! ...        steam or super critical region at this iteration (and was in
           ! ...        the two phase  region last time step) to a pressure and
           ! ...        enthalpy just outside the two phase/steam boundary.
           ! ...        Step along the line between the old P & H and the new
           ! ...        P & H until the point is just outside the boundary.
           ptemp = poldt(ic)
           htemp = holdt(ic)
           IF (ptemp >= 220.55e6_kdp) THEN
              hs = 20.86e9_kdp
           ELSE
              CALL phsbndy1(ptemp, dum1, hs)
           END IF
           aval = 0.0_kdp
           DO 
              IF (htemp > hs) EXIT
              aval = aval + ainc
              ptemp = aval*p(ic) + (1.0_kdp-aval)*poldt(ic)
              htemp = aval*h(ic) + (1.0_kdp-aval)*holdt(ic)
              IF (ptemp >= 220.55E6_kdp) THEN
                 hs = 20.86E9_kdp
              ELSE
                 CALL phsbndy1(ptemp, dum1, hs)
              END IF
           END DO
           p(ic) = ptemp
           h(ic) = htemp
        ELSEIF (indoldnr(ic) >= 3 .AND. ind(ic) == 2) THEN
           ! ... ****  PHASE CHANGE:   SUPERHEATED STEAM or SUPER CRITICAL to Two-PHASE ****
           ! ...      Adjust the pressure and enthalpy of any node that has gone from
           ! ...        superheated steam at the last N-R iteration to the two-phase region 
           ! ...        at this iteration
           ! ...        (and was in the superheated steam region last time step)
           ! ...        to a P and H just inside the two phase/steam boundary.
           ! ...        Step along the line between the old P & H and the new
           ! ...        P & H until the point is just inside the boundary.
           ptemp = poldt(ic)
           htemp = holdt(ic)
           IF (ptemp >= 220.55e6_kdp) THEN
              hs = 20.86e9_kdp
           ELSE
              CALL phsbndy1(ptemp, dum1, hs)
           END IF
           aval = 0.0_kdp
           DO 
              IF (htemp < hs .AND. ptemp < 220.55e6_kdp) EXIT
              aval = aval + ainc
              ptemp = aval*p(ic) + (1.0_kdp-aval)*poldt(ic)
              htemp = aval*h(ic) + (1.0_kdp-aval)*holdt(ic)
              IF (ptemp >= 220.55e6_kdp) THEN
                 hs = 20.86e9_kdp
              ELSE
                 CALL phsbndy1(ptemp, dum1, hs)
              END IF
           END DO
           p(ic) = ptemp
           h(ic) = htemp
        ELSE
           ! ...    write P & H of any node that has undergone a phase change
           ! ...         not covered by adjustments above
           IF (ioptpr(3) >= 1) WRITE(fuclog, 9010) '  P & H of node: ',  &
                p(ic), h(ic), 'no adjustment made'
           CYCLE
        END IF
        ! ...       write P & H of node after the adjustment for the phase change
        ! ...       also write the saturation H for water or steam at the new P
        IF (ioptpr(3) >= 1) THEN
           IF (p(ic) < 220.55e6_kdp) THEN
              CALL phsbndy1(p(ic),hw,hs)
              IF (indoldnr(ic) == 1 .OR. ind(ic) == 1) THEN
                 WRITE(fuclog, 9005) '               to ',p(ic),h(ic),  &
                      '; At this P, HW=', hw
              ELSE
                 WRITE(fuclog, 9005) '               to ',p(ic),h(ic),  &
                      '; At this P, HS=', hs
              END IF
           END IF
        END IF
        CYCLE
     END IF
  END DO
END SUBROUTINE phasechg
