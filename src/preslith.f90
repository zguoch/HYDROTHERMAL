SUBROUTINE preslith
  !     Purpose: To compute a lithostatic pressure gradient
  !     Method of calculating lithostatic pressure: Assign the "constant"
  !        pressure value input to the top cell in each column that is in
  !        the domain.  Working down the column of nodes, the pressure is
  !        estimated by from the saturated density of each node.  The saturated
  !        density is the density of rock * (1-porosity) * g * h + density of
  !        the liquid * porosity * g * h.
  !     NOTE: This subroutine does not account for 2-phase conditions in
  !        calculating the lithostatic pressure.
  USE machine_constants, ONLY: kdp
  USE f_units
!  USE bc
  USE i_c
  USE mesh
  USE parameters
  USE variables
  IMPLICIT NONE
  LOGICAL :: errflag
  INTEGER :: i, ic, icup, ik, j, k
  REAL(kind=kdp) :: db, dbu, dum1, dum2, dum3, dum4, dum5, dum6,  &
       dum7, dum8, dw, dwu
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5 $//$Date: 2002/12/19 19:53:53 $'
  !     ------------------------------------------------------------------
  !...
  DO  j = 1, ny
     DO  i = 1, nx
        icup = 0
        DO  k = nz, 1, -1
           ic = (k-1)*nx + i + (j-1)*nx*nz
           ik = (k-1)*nx + i
           !---skip past nodes not in the domain
           IF (np(ik,j) == 0) THEN
              plith(ic) = 0._kdp
              CYCLE
           END IF
           !---assign the fluid pressure in the top cell in each column
           !---to the lithostatic pressure
           IF (icup == 0) THEN
              plith(ic) = p(ic)
           ELSE
              ! ... use pressure-enthalpy function to calc fluid density in overlying block
              CALL tblookup(plith(icup), h(icup), dwu, dum1, dum2, dum3, dum4, dum5,  &
                   dum6, dum7, dum8, errflag)
              !---approx. litho pressure at node using water density of above node
              dbu = (1._kdp-phi(icup))*df(icup) + phi(icup)*dwu
              db = (1._kdp-phi(ic))*df(ic) + phi(ic)*dwu
              plith(ic) = plith(icup) + 0.5_kdp*dz(k+1)*grav*dbu + 0.5_kdp*dz(k)*grav*db
              ! ... recalculate density of water using new lithostatic pressure
              CALL tblookup(plith(ic), h(ic), dw, dum1, dum2, dum3, dum4, dum5, dum6,  &
                   dum7, dum8, errflag)
              ! ... re-approximate pressure at node
              db = (1._kdp-phi(ic))*df(ic) + phi(ic)*dw
              plith(ic) = plith(icup) + 0.5_kdp*dz(k+1)*grav*dbu + 0.5_kdp*dz(k)*grav*db
              ! ... once again recalculate water density for the given cell
              CALL tblookup(plith(ic), h(ic), dw, dum1, dum2, dum3, dum4, dum5, dum6,  &
                   dum7, dum8, errflag)
              ! ... calculate final pressure
              db = (1._kdp-phi(ic))*df(ic) + phi(ic)*dw
              plith(ic) = plith(icup) + 0.5_kdp*dz(k+1)*grav*dbu + 0.5_kdp*dz(k)*grav*db
              IF(errflag) THEN
                 WRITE (fustdout,9005)  &
                      '***** Pressure-Enthalpy table error for cell no. ',ic,' *****'
                 WRITE (fuclog,9005)  &
                      '***** Pressure-Enthalpy table error for cell no. ',ic,' *****'
9005             FORMAT(a,i6,a)
              END IF
           END IF
           icup = ic
        END DO
     END DO
  END DO
END SUBROUTINE preslith
