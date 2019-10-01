SUBROUTINE presboil(p,h,tc,ii,jj)
  ! ... Purpose: To compute a boiling point temperature profile with depth 
  ! ...          as a function of pressure and enthalpy
  ! ...           (temperature) profile assuming no steam in the column
  ! ... Method of calculating hydrostatic pressure profile: Assign the input
  ! ...    pressure value to the top cell in each column that is in
  ! ...    the domain.  Estimate the next pressure below by adding 
  ! ...    density of water*g*delz to the pressure of above cell.  
  ! ...    Calculate density and enthalpy of water in cell below from 
  ! ...    the saturated water curve.  Correct the pressure in the cell below.
  ! ...    Iterate the correction five times.
  ! ... ***should iterate to convergence***
  USE machine_constants, ONLY: kdp
  USE f_units
  USE bc
  USE parameters       !** for flint
  USE control
  USE mesh
!!$  USE parameters
!!$  USE variables
  IMPLICIT NONE
  REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: p, h, tc
  INTEGER, INTENT(IN) :: ii, jj
  !
  LOGICAL :: errflag
  INTEGER :: i, i1, i2, ic, icup, ik, j, j1, j2, k, lll, mt
  REAL(KIND=kdp) :: dum0, dum1, dum2, dum3, dum4, dum5, dum6, dum7,  &
       dum8, dum9, duma, dumb, dumc, dumd, dwavg, dw, dwu, pavg, pdwp
  REAL(KIND=kdp), PARAMETER :: pcrit=220.55e6_kdp
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.9 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  ! ... check that ptop < pcrit
  IF (ptop >= pcrit) THEN
     WRITE(fustdout,9005) ptop, pcrit
     WRITE(fuclog,9005) ptop, pcrit
9005 FORMAT (' ***** ERROR: ',  &
          'Problem defining boiling point temperature with depth conditions *****' , /, 11X,  &
          'Input pressure along top greater than critical point pressure', /,  &
          11X, 'Ptop=', 1pg12.4, 4X, 'Pcritical=', 1pg12.4)
     ierr(73) = .TRUE.
     RETURN
  END IF
  ! ... determine which columns to define boiling point temperature profile
  IF (ii == 0 .AND. jj == 0) THEN          ! ... do all the columns
     i1 = 1
     i2 = nx
     j1 = 1
     j2 = ny
  ELSE                                     ! ... do the specified column      
     i1 = ii
     i2 = ii
     j1 = jj
     j2 = jj
  END IF
  DO  j=j1,j2
     DO  i=i1,i2
        icup = 0
        !$$        mt=(j-1)*Nx+i
        DO  k = nz,1,-1
           ic = (k-1)*nx + i + (j-1)*nx*nz
           ik = (k-1)*nx + i
           ! ... skip past nodes not in the domain
           IF (np(ik,j) == 0) THEN
              p(ic) = 0._kdp
              CYCLE
           END IF
!!$           !..      if(ptopx(mt).ge.0._kdp) then
!!$           !..C---check that Ptop < pcrit
!!$           !..      IF (Ptopx(mt).GE.pcrit) THEN
!!$           !..        WRITE(fustdout, 9005) Ptopx(mt), pcrit
!!$           !..        WRITE(fuclog, 9005) Ptopx(mt), pcrit
!!$           !..        STOP
!!$           !..c              check that PTOP is with limits of the table
!!$           !..             elseif (Ptopx(mt).GT.1.0D10 .OR. Ptopx(mt).LT.0.5D6) THEN
!!$           !..                  WRITE(fustdout, 9055) Ptopx(mt)
!!$           !..                  WRITE(fuclog, 9055) Ptopx(mt)
!!$           !.. 9055 FORMAT (/, ' ***** STOP - ',
!!$           !..     &        'Error defining boiling point with depth conditions *****'
!!$           !..     &        , /, 11X,
!!$           !..     &        'Input pressure along top outside of table range', /,
!!$           !..     &        11X, 'Ptop=', 1PG12.4)
!!$           !..                  GOTO 30
!!$           !..                ENDIF
           IF (icup == 0) THEN
!!$              !..              P(ic) = Ptopx(mt)
              p(ic) = ptop        ! ... assign input pressure to top cell in domain 
              CALL phsbndy1(p(ic), h(ic), dum1)      ! ... Calculate enthalpy for saturated water
           ELSE
              ! ... calculate water density in the overlying cell and
              ! ...      estimate water density for the cell face by estimating the
              ! ...      change in density w.r.t. depth from
              ! ...      0.5* delta z * dDEN/dP * DEN(upper block) * gravity
              ! ... use the saturated water curve to calculate density
              IF (p(icup) <= pcrit) THEN           ! ... only for subcritical pressure
                 CALL phsbndy3(p(icup), dum1, dum2, dum3, dum4,  &
                      dwu, pdwp, dum5, dum6, dum7, dum8, dum9, duma,  &
                      dumb, dumc, dumd)
                 dwavg = dwu + 0.25_kdp*(dz(k+1)+dz(k))*pdwp*dwu*grav  ! *** 0.5??
                 ! ... approximate the pressure at node using estimated water density for 
                 ! ...     cell face
                 p(ic) = p(icup) + dwavg*grav*0.5_kdp*(dz(k+1)+dz(k))
                 ! ... calculate water density for the cell face
                 ! ...      using the average of the first estimate of pressure and the 
                 ! ...      pressure for the upper cell
                 ! ... iterate 5 times
                 !****fix someday to tolerance test***
                 DO  lll = 1, 5
                    pavg = (p(icup)+p(ic))*0.5_kdp
                    IF (pavg >= pcrit) EXIT
                    CALL phsbndy3(pavg, dum1, dum2, dum3, dum4,  &
                         dwavg, dum0, dum5, dum6, dum7, dum8, dum9, duma,  &
                         dumb, dumc, dumd)
                    p(ic) = p(icup) + dwavg*grav*0.5_kdp*(dz(k+1)+dz(k))
                 END DO
              END IF
              IF (p(ic) >= pcrit) THEN
                 h(ic) = 22.0e9_kdp        ! ... set h to supercritical value
              ELSE                         ! ... calculate enthalpy of water for the given cell
                 !                               using the most recent estimate of pressure
                 CALL phsbndy1(p(ic), h(ic), dum1)
              END IF
           END IF
           ! ... Decrease the enthalpy slightly so that the node is in the
           ! ...      single phase water region rather than on the saturation curve
           h(ic) = h(ic) - 0.01e9_kdp
           ! ... Calculate temperature of water for the given cell
           CALL tblookup(p(ic), h(ic), dw, dum1, dum2, tc(ic), dum4, dum5, dum3,  &
                dum6, dum7, errflag)
           IF(errflag) THEN
              WRITE (fustdout,9015)  &
                   '***** ERROR: Pressure-Enthalpy out of range in table for cell no. ',ic,' *****'
              WRITE (fuclog,9015)  &
                   '***** ERROR: Pressure-Enthalpy out of range in table for cell no. ',ic,' *****'
9015          FORMAT(a,i6,a)
           END IF
           icup = ic
        END DO
     END DO
  END DO
END SUBROUTINE presboil
