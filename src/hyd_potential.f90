SUBROUTINE hyd_potential
  ! ... Purpose:  To compute the hydraulic potential (units of pressure)
  ! ...      for steam and water.  This is only done for printing/plotting purposes.
  ! ...      the potential is not used in any other part of the program.
  ! ...      Also computes the potentiometric head based on local fluid density
  ! ...      and taking atmospheric pressure as datum.
  USE machine_constants, ONLY: kdp
  USE mesh
  USE parameters, ONLY: grav, patm
  USE variables
  IMPLICIT NONE
  INTEGER :: i, ic, ik, j, k
  REAL(KIND=kdp) :: gravh=980.665_kdp
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3 $//$Date: 2007/04/10 23:14:34 $'
  !     ------------------------------------------------------------------
  ! ... Calculate potential and head for steam and water 
  ! ... This is not the usual hydrodynamic potential or head
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     ! ... if node is not in domain, set potentials to zero
     IF (np(ik,j) == 0) THEN
        wpot(ic) = 0._kdp
        spot(ic) = 0._kdp
        hwpot(ic) = 0._kdp
        hspot(ic) = 0._kdp
     ELSE
        IF (satnw(ic) > 0._kdp) THEN
           ! ... hydraulic potential and head of water
           wpot(ic) = p(ic) + denw(ic)*grav*zz(k)
           hwpot(ic) = (p(ic)-patm)/(denw(ic)*gravh) + zz(k)
        END IF
        IF (ind(ic) /= 5 .AND. satnw(ic) < 1._kdp) THEN
           ! ... hydraulic potential and head of steam, if not air-water cell
           spot(ic) = p(ic) + dens(ic)*grav*zz(k)
           hspot(ic) =  (p(ic)-patm)/(dens(ic)*gravh) + zz(k)
        END IF
     END IF
  END DO
END SUBROUTINE hyd_potential
