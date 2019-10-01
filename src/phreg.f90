SUBROUTINE phreg(p,h)
  ! ... PURPOSE:  To determine the phase region of a cell from its pressure and enthalpy:
  ! ...      water, 2-phase, superheated steam, supercritical, or air-water region
  ! ... Method: The enthalpy of saturated water and saturated steam are
  ! ...      calculated for the pressure of the cell.  The enthalpy
  ! ...      of the cell relative to these values determines the phase region.
  ! ...   Ind(IC): 1 - compressed water region
  ! ...            2 - two phase region
  ! ...            3 - superheated steam
  ! ...            4 - supercritical fluid (P>Pc and H>Hc)
  ! ...            5 - air-water unsaturated zone
  USE machine_constants, ONLY: kdp
  USE f_units
  USE bc
  USE parameters          !!$*** for flint
  USE control
  USE mesh
!!$  USE parameters
  USE variables, ONLY: ind, indoldnr
  IMPLICIT NONE
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
  !
  INTEGER :: i, ic, ik, j, k
  REAL(KIND=kdp) :: hs, hw, hx, px, rlx
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.12 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     ind(ic) = 0
     IF (np(ik,j) == 0) CYCLE     ! ... skip if node is not in domain
     ! ... Unconfined flow, air-water zone
     !***this test will not work for van Genuchten; not available
     IF(unconfined .AND. p(ic) < pwb) THEN
        ind(ic) = 5
        CYCLE
     END IF
     ! ... Confined flow or saturated unconfined flow
     ! ... Check that P & H fall within the limits of the table
     px = p(ic)
     hx = h(ic)
     IF (px < 0.5e6_kdp) THEN
        IF(unconfined) THEN
           itmcntl(20) = itmcntl(20)+1
           IF(MOD(itmcntl(20),1000) == 1) THEN
              WRITE (fustdout,9015)  &
                   '***** Pressure too low for limits of the table, for P= ',  &
                   px, ' (dyne/cm^2);','Low limit will be used; occurrence ',  &
                   itmcntl(20),' *****'
              WRITE (fuclog,9015)  &
                   '***** Pressure too low for limits of the table, for P= ',  &
                   px, ' (dyne/cm^2);','Low limit will be used; occurrence ',  &
                   itmcntl(20),' *****'
9015          FORMAT (a,1pg14.7,a/tr10,a,i8,a)
           END IF
        ELSE
           WRITE (fustdout,9015)  &
                ' *****ERROR: Pressure too low for limits of the table, for P= ',  &
                px, ' (dyne/cm^2) *****'
           WRITE (fuclog,9015)  &
                ' *****ERROR: Pressure too low for limits of the table, for P= ',  &
                px, ' (dyne/cm^2) *****'
           ierr(135) = .TRUE.
        END IF
     END IF
     IF (px > 1.0e10_kdp) THEN
        WRITE (fustdout, 9015)  &
             ' *****ERROR: Pressure too high for limits of the table, for P= ',px,  &
             '(dyne/cm^2) *****'
        WRITE (fuclog, 9015)  &
             ' *****ERROR: Pressure too high for limits of the table, for P= ',px,  &
             '(dyne/cm^2) *****'
        ierr(123) = .TRUE.
     END IF
     IF (hx > 52.0e9_kdp) THEN
        WRITE (fustdout, 9015)  &
             ' *****ERROR: Enthalpy too high for limits of the table, for H= ',hx,'(erg/g) *****'
        WRITE (fuclog, 9015)  &
             ' *****ERROR: Enthalpy too high for limits of the table, for H= ',hx,'(erg/g) *****'
        ierr(124) = .TRUE.
     END IF
     IF ((px < 240.e6_kdp .AND. hx < 0.01e9_kdp) .OR.  &
          (px >= 240.0e6_kdp .AND. hx < 0.5e9_kdp)) THEN
        WRITE (fustdout, 9006)  &
             ' *****ERROR: Enthalpy too low for limits of the table, for H= ' , hx, '(erg/g)',  &
             'and P= ',px,'(dyne/cm^2) *****'
        WRITE (fuclog, 9006)  &
             ' *****ERROR: Enthalpy too low for limits of the table, for H= ' , hx, '(erg/g)',  &
             'and P= ',px,'(dyne/cm^2) *****'
9006    FORMAT (a,1PG14.7,a/tr10,a,1pg14.7,a)
        ierr(132) = .TRUE.
     END IF
     IF (ierr(135) .OR. ierr(123) .OR. ierr(124) .OR. ierr(132)) THEN
        WRITE (fustdout, 9010)  &
             ' ***ERROR: P or H Outside Limits of Table '//  &
             ' at I, J, K ',i, j, k,' P= ',p(ic),' H= ',h(ic)
        WRITE (fuclog, 9010)  &
             ' ***ERROR: P or H Outside Limits of Table '//  &
             ' at I, J, K ',i, j, k,' P= ',p(ic),' H= ',h(ic)
9010    FORMAT (a,3I3/tr10,a,1PG14.5,a,1PG14.5)
        ierr(71) = .TRUE.
        CYCLE
     END IF
     ! ... check for super-critical pressure region
     IF (p(ic) >= 220.55e6_kdp) THEN
        IF (h(ic) > 20.86e9_kdp) THEN
           ind(ic) = 4
        ELSE
           ind(ic) = 1
        END IF
        CYCLE
     END IF
     ! ... quick check for pure steam or compressed water region
     IF (h(ic) > 28.04e9_kdp) THEN
        ind(ic) = 3
        CYCLE
     END IF
     IF (h(ic) < 3.0e9_kdp .OR. h(ic) < ((LOG10(p(ic))-5.617_kdp)*5.217e9_kdp)) THEN
        ind(ic) = 1
        CYCLE
     END IF
     ! ... else it is in the two-phase region
     ind(ic) = 2
     ! ... More accurate check for enthalpy location in phase diagram
     CALL phsbndy1(p(ic), hw, hs)     ! ... Calculate enthalpy of saturated water and steam
     ! ... For the first 3 NR-iterations, use the true phase boundaries
     IF (lnri <= 3) THEN
        IF (h(ic) < hw) ind(ic) = 1
        IF (h(ic) > hs) ind(ic) = 3
        CYCLE
     ELSE
        ! ... For the subsequent NR-iterations, relax the phase boundaries
        ! ...          by a small amount to help reduce needless phase changes
        rlx = 0.0001e9_kdp
        ! ...          for node that was in the water region last NR iteration
        IF (indoldnr(ic) == 1) THEN
           !***this should be a fraction not absolute amount
           IF (h(ic) < (hw+rlx)) THEN
              ind(ic) = 1
              IF (h(ic) >= hw .AND. ioptpr(3) >= 1) WRITE (fuclog, 9005)  &
                   ltime, lnri, '2-PHASE', 'WATER', i, j, k, p(ic), h(ic),  &
                   ' at P, HW=',hw
9005          FORMAT (i4, i3, 3X, 'Node just into ', a,  &
                   ' region, but still treated as ', a, ' at', 3I3, /, 12X,  &
                   'P & H of node = ', 1P, 2D13.6, 1X, a8, d13.6)

           END IF
           IF (h(ic) > hs) ind(ic) = 3
        END IF
        ! ...          for node that was in the 2-phase region last NR iteration
        IF (indoldnr(ic) == 2) THEN
           IF (h(ic) < (hw-rlx)) ind(ic) = 1
           IF (h(ic) >= (hw-rlx) .AND. h(ic) < hw .AND.  &
                ioptpr(3) >= 1) WRITE (fuclog, 9005) ltime, lnri,  &
                'WATER', '2-PHASE', i, j, k, p(ic), h(ic),  &
                ' at P, HW=', hw
           IF (h(ic) > (hs+rlx)) ind(ic) = 3
           IF (h(ic) <= (hs+rlx) .AND. h(ic) > hs .AND.  &
                ioptpr(3) >= 1) WRITE (fuclog, 9005) ltime, lnri,  &
                'STEAM', '2-PHASE', i, j, k, p(ic), h(ic),  &
                ' at P, HS=', hs
        END IF
        ! ...          for node that was in the steam region last NR iteration
        IF (indoldnr(ic) == 3) THEN
           IF (h(ic) < hw) ind(ic) = 1
           IF (h(ic) > (hs-rlx)) THEN
              ind(ic) = 3
              IF (h(ic) <= hs .AND. ioptpr(3) >= 1)  &
                   WRITE (fuclog, 9005) ltime, lnri, '2-PHASE', 'STEAM',  &
                   i, j, k, p(ic), h(ic), ' at P, HS=',hs
           END IF
        END IF
        ! ...          for node that was in the super critical region last NR iteration
        IF (indoldnr(ic) == 4) THEN
           IF (h(ic) < hw) ind(ic) = 1
           IF (h(ic) > hs) ind(ic) = 3
        END IF
     END IF
     IF (ind(ic) < 1 .OR. ind(ic) > 5) THEN     ! ... error; invalid index for active cell
        WRITE(fustdout,'(A,3I3,/,10X,A,I1)') ' ***ERROR: Invalid Phase Index for node '// 'I,J,K  ',  &
             i, j, k, ' IND= ',ind(ic)
        WRITE(fuclog,'(A,3I3,/,10X,A,I1)') ' ***ERROR: Invalid Phase Index for node '// 'I,J,K  ',  &
             i, j, k, ' IND= ',ind(ic)
        ierr(102) = .TRUE.
     END IF
  END DO
END SUBROUTINE phreg
