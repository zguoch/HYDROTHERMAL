SUBROUTINE pickup
  ! ... Purpose:  To write the time label and pressure and enthalpy
  ! ...    and temperature and porosity arrays of the last time step
  ! ...    completed
  ! ...    in the format that is used by the input file
  ! ...    so that the program can be re-started
  ! ...    Written to Out_restartdump.* file
  USE machine_constants, ONLY: kdp
  USE f_units
  USE mesh
  USE parameters, ONLY: beta, phi
  USE variables
  USE units
  IMPLICIT NONE
  CHARACTER(LEN=80), PARAMETER :: cfmt='(1P,10(g12.5)/(10(g12.5)))'
  INTEGER :: ic, icxf, icxl, j, k
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.6 $//$Date: 2003/10/17 20:41:20 $'
  !     ------------------------------------------------------------------
  !...
  IF (iyr == 2) THEN
     WRITE(fuorst,'(1P,E12.6,A)') year(time),' (yr)'
  ELSE
     WRITE(fuorst,'(1P,E12.6,A)') time,' (s)'
  END IF
  ! ... porosity
  IF (beta(1) > 0.0_kdp) THEN
     WRITE(fuorst,9015) 'Porosity (-)'
9015 FORMAT (a)
     WRITE(fuorst,9020)
9020 FORMAT ('FREE'/'TOP')
     DO  j = 1, ny
        DO  k = nz, 1, -1
           icxf = (k-1)*nx + 1 + (j-1)*nx*nz
           icxl = icxf + nx - 1
           WRITE(fuorst, cfmt) (phi(ic), ic=icxf, icxl)
        END DO
     END DO
  END IF
  ! ... pressure
  WRITE(fuorst,9015) 'Pressure (dyne/cm^2)'
  WRITE(fuorst,9020)
  DO  j = 1, ny
     DO  k = nz, 1, -1
        icxf = (k-1)*nx + 1 + (j-1)*nx*nz
        icxl = icxf + nx - 1
        WRITE(fuorst, cfmt) (p(ic), ic=icxf, icxl)
     END DO
  END DO
  ! ... enthalpy
  WRITE(fuorst,9015) 'Enthalpy (erg/g)'
  WRITE(fuorst,9020)
  DO  j = 1, ny
     DO  k = nz, 1, -1
        icxf = (k-1)*nx + 1 + (j-1)*nx*nz
        icxl = icxf + nx - 1
        WRITE(fuorst, cfmt) (h(ic), ic=icxf, icxl)
     END DO
  END DO
  ! ... temperature
  WRITE(fuorst,9015) 'Temperature (Deg.C)'
  WRITE(fuorst,9020)
  DO  j = 1, ny
     DO  k = nz, 1, -1
        icxf = (k-1)*nx + 1 + (j-1)*nx*nz
        icxl = icxf + nx - 1
        WRITE(fuorst,cfmt) (tc(ic), ic=icxf, icxl)
     END DO
  END DO
END SUBROUTINE pickup
