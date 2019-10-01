SUBROUTINE parmcalc(dumarray,kodparm)
  ! ... Purpose:  Compute the porous media properties for each rock type
  ! ...           where properties are defined as functions of temperature.
  USE machine_constants, ONLY: kdp
  USE mesh
  USE parameters
  USE variables, ONLY: tc
  IMPLICIT NONE
  REAL(kind=kdp), DIMENSION(nxyz), INTENT(IN OUT) :: dumarray
  INTEGER, INTENT(IN) :: kodparm
  ! ...   DUMARRAY - property array that is to be assigned values that are f(Temperature)
  ! ...   KODPARM  - Rock Property code
  ! ...          1 - Porosity
  ! ...          2 - X permeability
  ! ...          3 - Y permeability
  ! ...          4 - Z permeability
  ! ...          5 - Thermal Conductivity
  ! ...          6 - Specific Heat
  ! ...          7 - Rock Density
  ! ...          8 - Rock Compressibility
  !
  INTEGER :: ic, icdum, irx, lll, funct_no
  REAL(kind=kdp) :: slope
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5 $//$Date: 2002/03/18 23:28:28 $'
  !     ------------------------------------------------------------------
  loop30:  DO  icdum = 1, npiccons
     ic = npic(icdum)
     irx = icrxtype(ic)
     funct_no = irxftopt(kodparm,irx,1)
     IF (funct_no == 1) THEN      ! ... function 1 - parameter NOT a function of T
        dumarray(ic) = rxftparm(kodparm,irx,1)
     ELSE IF (funct_no == 2) THEN  ! ... function 2 - step function of T
        IF (tc(ic) < rxftparm(kodparm,irx,2)) THEN
           dumarray(ic) = rxftparm(kodparm,irx,1)
        ELSE
           DO  lll=1,irxftopt(kodparm,irx,2)-1
              IF (tc(ic) >= rxftparm(kodparm,irx,2*lll)) dumarray(ic) =  &
                   rxftparm(kodparm,irx,2*lll+1)
           END DO
        END IF
     ELSE IF (funct_no == 3 .OR. funct_no == 4) THEN 
        ! ... function 3 - log parameter is a linear function of T
        ! ...     or function 4 - parameter is a linear function of T
        IF (tc(ic) < rxftparm(kodparm,irx,2)) THEN    ! ...  T less than first input Temp point
           dumarray(ic) = rxftparm(kodparm,irx,1)
        ELSE IF (tc(ic) >= rxftparm(kodparm,irx,2*irxftopt(kodparm,irx,2))) THEN
           ! ...  T greater than last input Temp point
           dumarray(ic) = rxftparm(kodparm,irx,2*irxftopt(kodparm,irx,2 )-1)
        ELSE
           DO  lll=1,irxftopt(kodparm,irx,2)-1
              IF (tc(ic) >= rxftparm(kodparm,irx,2*lll) .AND. tc(ic)  &
                   < rxftparm(kodparm,irx,2*lll+2)) THEN
                 IF (funct_no == 3) THEN
                    slope = (LOG10(rxftparm(kodparm,irx,2*lll+1))  &
                         -LOG10(rxftparm(kodparm,irx,2*lll-1)))  &
                         /(rxftparm(kodparm,irx,2*lll+2) -rxftparm(kodparm,irx,2*lll))
                    dumarray(ic) = 10._kdp**(LOG10(rxftparm(kodparm,irx,2* lll-1))  &
                         +slope*(tc(ic)-rxftparm(kodparm,irx, 2*lll)))
                 ELSEIF (funct_no == 4) THEN
                    slope = (rxftparm(kodparm,irx,2*lll+1)  &
                         -rxftparm(kodparm,irx,2*lll-1))/(rxftparm(kodparm,irx,2*lll+2)  &
                         -rxftparm(kodparm,irx,2*lll))
                    dumarray(ic) = rxftparm(kodparm,irx,2*lll-1) + slope*(tc(ic)  &
                         -rxftparm(kodparm,irx,2*lll))
                 END IF
                 CYCLE loop30
              END IF
           END DO
        END IF
     END IF
  END DO loop30
END SUBROUTINE parmcalc
