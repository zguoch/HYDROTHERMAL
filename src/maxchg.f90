SUBROUTINE maxchg(wschgocr, pchgocr, hchgocr, iwschg, jwschg,  &
     kwschg, ipchg, jpchg, kpchg, ihchg, jhchg, khchg, atmcut)
  ! ... Purpose:  To determine the maximum change in water saturation and
  ! ...      the max % change in P and H over the current time step.  The
  ! ...      change in water saturation is computed only for pressures < 200 bar
  ! ... This is used for time step control
  USE machine_constants, ONLY: kdp
  USE f_units
  USE bc
  USE control
  USE mesh
  USE variables
  IMPLICIT NONE
  REAL(KIND=kdp), INTENT(OUT) :: wschgocr, pchgocr, hchgocr
  INTEGER, INTENT(OUT) :: iwschg, jwschg, kwschg, ipchg, jpchg, kpchg, ihchg, jhchg, khchg
  CHARACTER (LEN=1), INTENT(INOUT) :: atmcut     ! ... abbreviation of reason for cutting time step
  !
  INTEGER :: i, ic, ik, j, k
  REAL(KIND=kdp) :: awschgocr, apchgocr, ahchgocr, delth, deltp, deltsw
  REAL(KIND=kdp), PARAMETER :: faclimit=1.1_kdp
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.11 $//$Date: 2007/10/16 21:17:22 $'
  !     ------------------------------------------------------------------
  wschgocr = 0.0_kdp
  pchgocr = 0.0_kdp
  hchgocr = 0.0_kdp
  iwschg = 0
  jwschg = 0
  kwschg = 0
  ipchg = 0
  jpchg = 0
  kpchg = 0
  ihchg = 0
  jhchg = 0
  khchg = 0
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     ! ... Skip if node is a constant P and H cell or is not in active domain
     ! ...    or is a seeping cell
     IF(ibc(ik,j) == -1 .OR. ibc(ik,j) == 11 .OR.  &
          ((ibc(ik,j) == 30 .OR. ibc(ik,j) == 50) .AND. seeping(ik,j))) CYCLE
     ! ... Change in water saturation
     ! ... do not check saturation change when pressure > 200.0 bar
     IF (poldt(ic) <= 200.e6_kdp .AND. p(ic) <= 200.e6_kdp) THEN
        IF (ABS(satnw(ic)-swoldt(ic)) >= ABS(wschgocr)) THEN
           wschgocr = satnw(ic)-swoldt(ic)
           iwschg = i
           jwschg = j
           kwschg = k
        END IF
     END IF
     ! ... Change in pressure
     IF(unconfined) THEN           ! ... change in pressure, unscaled
        IF(ABS(p(ic)-poldt(ic)) >= ABS(pchgocr)) THEN
           pchgocr = p(ic)-poldt(ic)
           ipchg = i
           jpchg = j
           kpchg = k
        END IF
     ELSE                 ! ...  percent change in pressure (confined flow only)
        IF (ABS(100._kdp*(p(ic)-poldt(ic))/poldt(ic)) >= ABS(pchgocr)) THEN
           pchgocr = 100._kdp*(p(ic)-poldt(ic))/poldt(ic)
           ipchg = i
           jpchg = j
           kpchg = k
        END IF
     END IF
     ! ...  Percent change in enthalpy
     IF (ABS(100._kdp*(h(ic)-holdt(ic))/holdt(ic)) >= ABS(hchgocr)) THEN
        hchgocr = 100._kdp*(h(ic)-holdt(ic))/holdt(ic)
        ihchg = i
        jhchg = j
        khchg = k
     END IF
  END DO
  ! ... make sure absolute values for change do not equal 0
  awschgocr = MAX(ABS(wschgocr),1.e-25_kdp)
  apchgocr = MAX(ABS(pchgocr),1.e-25_kdp)
  ahchgocr = MAX(ABS(hchgocr),1.e-25_kdp)
  ! ... Check if any changes exceed the maximum allowed and
  ! ...     adjust the time step if necessary
  ! ... Time step is cut if max dependent variable changes are greater than
  ! ...      faclimit times specified limits
  ! ... ***Presently enforces 1.1 times the limit for cut in time step
  ! ... Check the change in water saturation
  deltsw = delt
  IF (awschgocr > faclimit*wschgmx) THEN
     deltsw = delt*wschgmx/awschgocr
     IF (ioptpr(3) >= 1) THEN
        WRITE (fustdout,9005) wschgocr, wschgmx, iwschg, jwschg, kwschg
        WRITE (fuclog,9005) wschgocr, wschgmx, iwschg, jwschg, kwschg
9005    FORMAT (tr5,'Time-step Saturation change of',f7.4,  &
             ' exceeds limit (',f6.4,'); I,J,K: ',3I3)
     END IF
  END IF
  ! ... Check the change in pressure (percent or absolute)
  deltp = delt
  IF (apchgocr > faclimit*pchgmxts) THEN
     deltp = delt*pchgmxts/apchgocr
     IF (ioptpr(3) >= 1) THEN
        IF(unconfined) THEN
           WRITE (fustdout,9110) 'Pressure', pchgocr, ipchg, jpchg, kpchg
           WRITE (fuclog,9110) 'Pressure', pchgocr, ipchg, jpchg, kpchg
9110       FORMAT(tr5,'Time-step ',a,' change of',1pg10.3,  &
                ' (dy/cm^2) exceeds limit; I,J,K: ',3I3)
        ELSE
           WRITE (fustdout, 9010) 'Pressure', pchgocr, ipchg, jpchg, kpchg
           WRITE (fuclog, 9010) 'Pressure', pchgocr, ipchg, jpchg, kpchg
9010       FORMAT(tr5,'Time-step ',a,'  % change of',f10.3,  &
                ' exceeds limit; I,J,K: ',3I3)
        END IF
     END IF
  END IF
  ! ... Check the change in enthalpy (percent)
  delth = delt
  IF (ahchgocr > faclimit*hchgmxts) THEN
     delth = delt*hchgmxts/ahchgocr
     IF (ioptpr(3) >= 1) THEN
        WRITE (fustdout,9010) 'Enthalpy', hchgocr, ihchg, jhchg, khchg
        WRITE (fuclog,9010) 'Enthalpy', hchgocr, ihchg, jhchg, khchg
     END IF
  END IF
  IF (deltsw < delt .OR. deltp < delt .OR. delth < delt) THEN
     ! ... The maximum change in one or more of water saturation, P, or H exceeded
     ! ...     the specified limits.
     ! ... use the minimum reduced time step
     delt = MIN(deltsw, deltp, delth)
     IF (delt == deltsw) THEN
        itmcntl(11) = itmcntl(11) + 1
        atmcut = 'W'
     ELSEIF (delt == deltp) THEN
        itmcntl(12) = itmcntl(12) + 1
        atmcut = 'P'
     ELSEIF (delt == delth) THEN
        itmcntl(13) = itmcntl(13) + 1
        atmcut = 'H'
     END IF
     tmcut = .TRUE.
  END IF
END SUBROUTINE maxchg
