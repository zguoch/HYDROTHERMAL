SUBROUTINE press(ndenw,h,tc,p)
  ! ... Purpose: To compute a hydrostatic pressure distribution with depth for liquid
  ! ...          water only.
  ! ... Method of calculating hydrostatic pressure profile: Assign the input
  ! ...    pressure value to the top cell in each column that is in
  ! ...    the domain.  Estimate the next pressure below by adding 
  ! ...    density of water*g*delz to the pressure of above cell.  
  ! ...    Calculate density from pressure and enthalpy (temperature)
  ! ...    Correct the pressure in the cell below.
  ! ...    Iterate the correction two times.
  ! ... ***should iterate to convergence***
  USE machine_constants, ONLY: kdp
  USE f_units
  USE bc
  USE mesh
  USE parameters, ONLY: grav
  USE steam_table_funct
  USE variables, ONLY: ktop
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ndenw
  ! ...   NDENW: indicates method for calculating water density
  ! ...         1 - use temperature and pressure in subroutine TEMPRESS
  ! ...         2 - use enthalpy and pressure in subroutine TBLOOKUP
  ! ...         3 - mixed temperature and enthalpy data. use appropriate routine
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: h, tc, p
  !
  LOGICAL :: errflag
  INTEGER :: i, ic, ickp, ij, ik, j, k
  REAL(KIND=kdp) :: dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dw, dwu
  INTERFACE
     SUBROUTINE tblookup(px, hx, da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp, errflag)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), INTENT(IN) :: hx, px
       REAL(KIND=kdp), INTENT(OUT) :: da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp
       LOGICAL, INTENT(INOUT) :: errflag
     END SUBROUTINE tblookup
  END INTERFACE
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.10 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  !...
  DO  j = 1, ny
     DO  i = 1, nx
        ij = (j-1)*nx + i
        DO  k = nz,1,-1
           ic = (k-1)*nx + i + (j-1)*nxz
           ickp = k*nx + i + (j-1)*nxz
           ik = (k-1)*nx + i
           errflag = .FALSE.
           IF (k > ktop(ij)) THEN
              p(ic) = 0._kdp              ! ... Null value for nodes not in the domain
           ELSEIF(k == ktop(ij)) THEN
              ! ... Assign the input pressure to the uppermost active cell in each column
              p(ic) = ptopa(ij)
           ELSE
              ! ... Calculate water density in the current overlying cell and
              ! ...      estimate water density for the cell in question
              IF (ndenw == 1) THEN
                 CALL tempress(tc(ickp), p(ickp), dwu, dum1)     ! ... calculate density
                 p(ic) = p(ickp) + 0.5_kdp*grav*dwu*(dz(k+1) + dz(k))    ! ... Estimate pressure
                 ! ...                     at current node using density of above node
                 CALL tempress(tc(ic), p(ic), dw, dum1)    ! ... Calculate density from estimated 
                 ! ...                                          pressure
              ELSEIF (ndenw == 2) THEN
                 CALL tblookup(p(ickp), h(ickp), dwu, dum1, dum2, dum3, dum4, dum5, dum6,  &
                      dum7, dum8, errflag)          ! ... calculate density
                 p(ic) = p(ickp) + 0.5_kdp*grav*dwu*(dz(k+1) + dz(k))    ! ... Estimate pressure
                 ! ...                     at current node using density of above node
                 CALL tblookup(p(ic), h(ic), dw, dum1, dum2, dum3, dum4, dum5, dum6,  &
                      dum7, dum8, errflag)         ! ... Calculate density from estimated pressure
              ELSEIF (ndenw == 3) THEN
                 IF (tc(ickp) <= -99._kdp) THEN     ! ... assume that enthalpy value input
                    CALL tblookup(p(ickp), h(ickp), dwu, dum1, dum2, dum3, dum4, dum5,  &
                         dum6, dum7, dum8, errflag)
                 ELSE                               ! ... assume that temperature value input
                    CALL tempress(tc(ickp), p(ickp), dwu, dum1)
                 END IF
                 p(ic) = p(ickp) + 0.5_kdp*grav*dwu*(dz(k+1) + dz(k))    ! ... Estimate pressure
                 ! ...                     at current node using density of above node
                 IF (tc(ic) <= -99._kdp) THEN     ! ...  assume that enthalpy value input
                    CALL tblookup(p(ic), h(ic), dw, dum1, dum2, dum3, dum4, dum5, dum6,  &
                         dum7, dum8, errflag)
                 ELSE                               ! ... assume that temperature value input
                    CALL tempress(tc(ic), p(ic), dw, dum1)
                 END IF
              END IF
              ! ... Recalculate pressure at node using the corrected water density for the node
              p(ic) = p(ickp) + 0.5_kdp*grav*(dwu*dz(k+1) + dw*dz(k))
              ! ... Recalculate water density for the current cell
              IF(ndenw == 1) THEN
                 CALL tempress(tc(ic), p(ic), dw, dum1)
              ELSEIF(ndenw == 2) THEN
                 CALL tblookup(p(ic), h(ic), dw, dum1, dum2, dum3, dum4, dum5, dum6,  &
                      dum7, dum8, errflag)
              ELSEIF(ndenw == 3) THEN
                 IF (tc(ic) <= -99._kdp) THEN     ! ...  assume that enthalpy value input
                    CALL tblookup(p(ic), h(ic), dw, dum1, dum2, dum3, dum4, dum5, dum6,  &
                         dum7, dum8, errflag)
                 ELSE                               ! ... assume that temperature value input
                    CALL tempress(tc(ic), p(ic), dw, dum1)
                 END IF
              END IF
              ! ... Recalculate pressure at current node
              p(ic) = p(ickp) + 0.5_kdp*grav*(dwu*dz(k+1) + dw*dz(k))
           END IF
           IF(errflag) THEN
              WRITE (fustdout,9005)  &
                   '***** ERROR: Pressure-Enthalpy out of range for table for cell no. ',ic,' or',ickp,' *****'
              WRITE (fuclog,9005)  &
                   '***** ERROR: Pressure-Enthalpy out of range for table for cell no. ',ic,' or',ickp,' *****'
9005          FORMAT(a,i6,a,i6,a)
           END IF
        END DO
     END DO
  END DO
END SUBROUTINE press
