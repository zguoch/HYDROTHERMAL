SUBROUTINE enthrock
  ! ... Purpose:  To compute the enthalpy of the rock matrix and the
  ! ...     derivative of the rock-matrix heat capacity
  USE machine_constants, ONLY: kdp
  USE mesh
  USE parameters
  USE variables
  IMPLICIT NONE
  INTEGER :: ic, icdum, irx, irxok, lll
  REAL(kind=kdp) :: bintrcpt, slope
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5 $//$Date: 2006/09/06 21:04:28 $'
  !     ------------------------------------------------------------------
  !...
  ! ... determine the rock type number for the first active node
  irxok = icrxtype(npic(1))
  loop30:  DO  icdum = 1, npiccons
     ic = npic(icdum)
     irx = icrxtype(ic)
     ! ... option 1 - Calculate enthalpy of rock for heat capacity NOT a function of T
     ! ...     that is where:
     ! ...       rock type = 0 (ie rock type input not used), or
     ! ...       Cp input as a constant, or
     ! ...       Cp is NOT a function of T (option 1 for IRXFTOPT) for this rock type
     IF (irxok == 0 .OR. irxftopt(6,irxok,1) == 0 .OR. irxftopt(6,irx,1) == 1) THEN
        hrock(ic) = phfwt(ic)*(tc(ic)+273.15_kdp)
        phfwtdt(ic) = 0.0_kdp
        CYCLE loop30
     END IF
     ! ... option 2 - Cp is a step function of T
     IF (irxftopt(6,irx,1) == 2) THEN
        phfwtdt(ic) = 0.0_kdp
        hrock(ic) = 0.0_kdp
        ! ... if T<first Temperature step
        IF (tc(ic) < rxftparm(6,irx,2)) THEN
           hrock(ic) = rxftparm(6,irx,1)*(tc(ic)+273.15_kdp)
           phfwtdt(ic) = 0.0_kdp
           CYCLE loop30
        END IF
        hrock(ic) = rxftparm(6,irx,1)*(rxftparm(6,irx,2)+273.15_kdp)
        ! ... do all but the first and last temperature steps
        IF (irxftopt(6,irx,2) >= 3) THEN
           DO  lll = 2, irxftopt(6,irx,2) - 1
              IF (tc(ic) < rxftparm(6,irx,2*lll)) THEN
                 hrock(ic) = hrock(ic) + rxftparm(6,irx,2*lll-1)*(tc(ic)-rxftparm(6,irx,2*lll-2))
                 CYCLE loop30
              END IF
              hrock(ic) = hrock(ic) + rxftparm(6,irx,2*lll-1)*  &
                   (rxftparm(6,irx,2*lll)-rxftparm(6,irx,2*lll-2))
           END DO
        END IF

        IF (tc(ic) >= rxftparm(6,irx,2*irxftopt(6,irx,2)-2)) THEN
           ! ...  Temperature >= last step
           hrock(ic) = hrock(ic) + rxftparm(6,irx,2*irxftopt(6,irx,2)-1)*  &
                (tc(ic)-rxftparm(6,irx,2*irxftopt(6,irx,2)-2))
           CYCLE loop30
        END IF
     END IF
!!$     ! ... option 3 - log Cp is a linear function of T
!!$     ! ...   ??? Note: this option not yet working for Cp
!!$        IF (Irxftopt(6,irx,1).EQ.3 .OR.
!!$              IF (Irxftopt(6,irx,1).EQ.3) THEN
!!$                slope = (LOG10(Rxftparm(6,irx,2*lll+1)) -
!!$     &                    LOG10(Rxftparm(6,irx,2*lll-1)))
!!$     &              / (Rxftparm(6,irx,2*lll+2) -
!!$     &                  Rxftparm(6,irx,2*lll))
!!$                DUMARRAY(ic) =
!!$     &            10.0_kdp**(LOG10(Rxftparm(6,irx,2*lll-1)) +
!!$     &              slope*(TC(IC)-Rxftparm(6,irx,2*lll)))
!!$              ENDIF
     ! ... option 4 - Cp is a linear function of T
     IF (irxftopt(6,irx,1) == 4) THEN
        hrock(ic) = 0.0_kdp
        ! ... if T<first input Temp point
        IF (tc(ic) < rxftparm(6,irx,2)) THEN
           hrock(ic) = rxftparm(6,irx,1)*(tc(ic)+273.15_kdp)
           phfwtdt(ic) = 0.0_kdp
           CYCLE loop30
        END IF
        hrock(ic) = rxftparm(6,irx,1)*(rxftparm(6,irx,2)+273.15_kdp)
        ! ... check for temperature between the first and last points
        DO  lll = 1, irxftopt(6,irx,2) - 1
           slope = (rxftparm(6,irx,2*lll+1)- rxftparm(6,irx,2*lll-1))/  &
                (rxftparm(6,irx,2*lll+2)- rxftparm(6,irx,2*lll))
           bintrcpt = rxftparm(6,irx,2*lll-1) - slope*rxftparm(6,irx,2*lll)
           IF (tc(ic) < rxftparm(6,irx,2*lll+2)) THEN
              hrock(ic) = hrock(ic) + slope/2.0_kdp*((tc(ic)+273.15_kdp)**2 -  &
                   (rxftparm(6,irx,2*lll)+273.15_kdp)**2) + bintrcpt*  &
                   (tc(ic)-rxftparm(6,irx,2*lll))
              phfwtdt(ic) = slope
              CYCLE loop30
           END IF
           hrock(ic) = hrock(ic) + slope/2.0_kdp*  &
                ((rxftparm(6,irx,2*lll+2)+273.15_kdp)**2 -  &
                (rxftparm(6,irx,2*lll)+273.15_kdp)**2) + bintrcpt*  &
                (rxftparm(6,irx,2*lll+2)-rxftparm(6,irx,2*lll))
        END DO
        IF (tc(ic) >= rxftparm(6,irx,2*irxftopt(6,irx,2))) THEN
           ! ...  T >= last input Temperature point
           hrock(ic) = hrock(ic) + rxftparm(6,irx,2*irxftopt(6,irx,2)-1)*  &
                (tc(ic)-rxftparm(6,irx,2*irxftopt(6,irx,2)))
           phfwtdt(ic) = 0.0_kdp
        END IF
     END IF
  END DO loop30
END SUBROUTINE enthrock

