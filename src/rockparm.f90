SUBROUTINE rockparm
  ! ... Purpose:  To compute the permeability, porosity, thermal conductivity,
  ! ...           heat capacity, density, compressibility of porous matrix
  ! ...           for selected rock types as a function of temperature.
  !
  USE machine_constants, ONLY: kdp
  USE mesh
  USE parameters
  IMPLICIT NONE
  INTEGER :: ic, icdum, irx, irxok
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5 $//$Date: 2003/10/15 22:38:56 $'
  !     ------------------------------------------------------------------
  ! ...    determine the rock type number for the first active node
  irxok = icrxtype(npic(1))
  ! ... Porosity
  IF (irxftopt(1,irxok,1) > 0) CALL parmcalc(phi,1)
  ! ... X Permeability
  IF (irxftopt(2,irxok,1) > 0) CALL parmcalc(xk,2)
  ! ... Y Permeability
  IF (irxftopt(3,irxok,1) > 0) THEN
     CALL parmcalc(yk,3)
  ELSE IF (irxftopt(2,irxok,1) > 0 .AND. ykrxkxfac(1) /= 0._kdp) THEN
     DO icdum = 1,npiccons
        ic = npic(icdum)
        irx = icrxtype(ic)
        yk(ic) = ykrxkxfac(irx)*xk(ic)
     END DO
  END IF
  ! ... Z Permeability
  IF (irxftopt(4,irxok,1) > 0) THEN
     CALL parmcalc(zk,4)
  ELSE IF (irxftopt(2,irxok,1) > 0 .AND. zkrxkxfac(1) /= 0._kdp) THEN
     DO icdum = 1,npiccons
        ic = npic(icdum)
        irx = icrxtype(ic)
        zk(ic) = zkrxkxfac(irx)*xk(ic)
     END DO
  END IF
  ! ... Thermal Conductivity
  IF (irxftopt(5,irxok,1) > 0) CALL parmcalc(xkc,5)
  ! ... Specific Heat (Heat Capacity)
!!$  ! ...       calculate derivate of Heat Capacity w.r.t. Temperature
!!$  ! ...       if first time step then set to 0  *** this is disabled
!!$  ! ...   IF (Ltime.GT.1) CALL PARMCAL2 (Phfwt,Phfwtoldt,Phfwtdt,6)
  IF (irxftopt(6,irxok,1) > 0) CALL parmcalc(phfwt,6)
  ! ... Rock Density
  IF (irxftopt(7,irxok,1) > 0) CALL parmcalc(df,7)
  ! ... Rock Compressibility 
  IF (irxftopt(8,irxok,1) > 0) CALL parmcalc(beta,8)
END SUBROUTINE rockparm
