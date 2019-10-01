SUBROUTINE enthtemp(p,tc,h)
  ! ... Purpose: To compute enthalpy from temperature and pressure
  ! ...      Used when temperature was input or when enthalpy is to be set to the 
  ! ...      saturated water value for the pressure
  ! ... Method of calculation: Call subroutine tempress
  ! ...    which uses the steam tables to calculate enthalpy from
  ! ...    temperature and pressure.
  ! ...    N-R iteration:  new enthalpy = old enthalpy - delta enthalpy
  ! ...       where:  delta enth = delta T / slope of function at old enth
  ! ...    If pressure was read as a constant, the hydrostatic profile is
  ! ...       first calculated using the temperature data.  The hydrostatic
  ! ...       pressure is later recalculated using the computed enthalpy.
  USE machine_constants, ONLY: kdp
  USE control
  USE mesh
  USE steam_table_funct
!!$  USE variables, ONLY: h, p
  IMPLICIT NONE
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: p, tc, h
  !
  INTEGER :: i, ic, ik, j, k
  REAL(KIND=kdp) :: dum1, hh
  INTERFACE
     SUBROUTINE press(ndenw,h,tc,p)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: ndenw
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: h, tc, p
     END SUBROUTINE press
  END INTERFACE
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.6 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  IF (nz > 1 .AND. kod(11) == 60 .AND. kodt > 0 .AND. kod(10)  < 4) THEN
     ! ...  a hydrostatic profile is to be calculated for pressure,
     ! ...      and enthalpy is not from saturated water curve.
     ! ...  compute hydrostatic profile using temperature data
     SELECT CASE (kodt)
     CASE (1)
        CALL press(1,h,tc,p)     ! ...  temperature read for all nodes
     CASE (2)
        CALL press(3,h,tc,p)     ! ...  temperature read for a subset of nodes
     END SELECT
  END IF
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF (np(ik,j) == 0) THEN           ! ... node not in domain; load 0; skip to end of loop
        h(ic) = 0._kdp
        CYCLE
     END IF
     IF((kodt == 2 .OR. kod(10) == 77) .AND. tc(ic) <= -99._kdp) CYCLE
     ! ... assume that enthalpy value was input for that node
     ! ...      instead of temperature so skip enthalpy calculation
     IF (kod(10) == 77 .AND. h(ic) == -77._kdp) THEN
        CALL phsbndy1(p(ic), h(ic), dum1)      ! ... calculate enthalpy for saturated water
     ELSE
        CALL tempress(tc(ic), p(ic), dum1, h(ic))  ! ... calculate enthalpy from temperature
     END IF
  END DO
END SUBROUTINE enthtemp
