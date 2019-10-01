SUBROUTINE init_wr_bc(initial)
  ! ... Purpose: To set up and write initial or new values for the
  ! ...    boundary conditions:
  ! ...           heat flux at base
  ! ...           precipitation flux at top surface
  USE machine_constants, ONLY: kdp
  USE math_constants
  USE f_units
  USE control, ONLY: ierr, ioptpr, kod
  USE bc
  USE mesh
  USE parameters
  USE steam_table_funct
  USE units, ONLY: unitfac, unitlabel
  USE variables, ONLY: ktop, rsq, visw
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: initial
  ! ...   Initial - .TRUE.= initial values input
  ! ...             .FALSE.= new values input
  !
  CHARACTER(LEN=80) :: cfmt
  CHARACTER(LEN=44), PARAMETER :: cfmt1="('Row  1 -> ',1P,6(g11.4)/(10X,6(g11.4))) ",  &
       cfmt2="('Row  1 -> ',1P,11(g11.4)/(10X,11(g11.4)))",  &
       cfmt3="('Row  1 -> ',1P,13(g11.4)/(10X,13(g11.4)))"
  INTEGER :: i, icpt, ij, ijxf, ijxl, ik, ikpt, j, k
  ! ... ------------------------------------------------------------------
  ! ...   Kod() =1 - array input as variable value array
  ! ...         =2 - array input as constant value array
  ! ...         =3 - array input as values constant for each row
  ! ...         =4 - array input as values constant for each column
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.7 $//$Date: 2007/08/03 23:42:47 $'
  ! ... ------------------------------------------------------------------
  !...
  ! ... assign format for writing fluxes depending on chosen output width
  IF (ioptpr(2) == 0) THEN
     cfmt = cfmt1
  ELSE IF (ioptpr(2) == 1) THEN
     cfmt = cfmt2
  ELSE
     cfmt = cfmt3
  END IF
  ! ... Basal Heat Flux (Conduction)
  ! ... KOD(9): 1 for spatially variable flux; 2 for uniform flux
  nhcond = 0
  IF (kod(9) > 0) THEN
     ! ... Set flags in IBC array for all the active heat conduction nodes
     !...****only works for flat bottom regions at k=1. Fix later, someday.
     DO  j=1,ny
        k = 1          ! .. only at bottom of global mesh
        DO  i=1,nx
           ij = (j-1)*nx + i
           ik = (k-1)*nx + i
           IF(ABS(cdtn(ij)) > 0._kdp ) THEN
              IF(ibc(ik,j) == 0) THEN
                 ibc(ik,j) = 02
                 nhcond = nhcond + 1
              END IF
           END IF
        END DO
     END DO
     IF (initial) THEN
        WRITE(fupdef,9005) 'Basal Heat Flux', unitlabel(9)
     ELSE
        WRITE(fupdef,9005) 'New values for Basal Heat Flux', unitlabel(9)
     END IF
     ! ... convert to user output units and write array
     IF (kod(9) == 2) THEN    ! ... uniform flux
        WRITE(fupdef,9025) cdtn(1)/unitfac(9)
     ELSE                     ! ... spatial distribution
        DO  j = 1,ny
           WRITE(fupdef, '(tr1)')
           WRITE(fupdef,9030) '--- Slice',j,' ---'
           ijxf = (j-1)*nx + 1
           ijxl = j*nx
           WRITE(fupdef,cfmt) (cdtn(ij)/unitfac(9), ij=ijxf,ijxl)
        END DO
     END IF
     IF (irad) THEN     ! ... basal surface area for cylindrical coordinates
        WRITE(fupdef, '(/tr10,A/tr15,A)')  &
             'Basal Surface Area per Cell for Cylindrical Coordinates (cm^2)',  &
             ' (Calculated from squared nodal cell boundaries)'
        WRITE(fupdef, '(1X)')
        WRITE(fupdef, cfmt) (pi*(rsq(i+1)-rsq(i)), i=1,nx)
     ELSE        ! ... basal surface area for Cartesian coordinates
        WRITE(fupdef, '(/tr10,A)')  &
             'Basal Surface Area per Cell for Cartesian Coordinates (cm^2)'
        IF (ABS(kod(12)) == 2 .AND. ABS(kod(13)) == 2) THEN
           WRITE(fupdef,9025) dx(1)*dy(1)
        ELSE
           DO  j=1,ny
              WRITE(fupdef, '(1X)')
              WRITE(fupdef,9030) '--- Slice',j,' ---'
              WRITE(fupdef,cfmt) (dx(i)*dy(j), i=1,nx)
           END DO
        END IF
     END IF
     ! ... basal heat input per node
     WRITE(fupdef, '(/tr10,A)') 'Heat Input flow rate per cell (erg/s)'
     IF (irad) THEN
        WRITE(fupdef,'(tr1)')
        WRITE(fupdef,cfmt) (cdtn(i)*pi*(rsq(i+1)-rsq(i)), i=1,nx)
     ELSE IF (kod(9) == 2 .AND. ABS(kod(12)) == 2 .AND. ABS(kod(13)) == 2) THEN
        WRITE(fupdef,9025) cdtn(1)*dx(1)*dy(1)
     ELSE
        DO  j=1,ny
           WRITE(fupdef,'(tr1)')
           WRITE(fupdef,9030) '--- Slice',j,' ---'
           WRITE(fupdef,cfmt) (cdtn((j-1)*nx+i)*dx(i)*dy(j), i=1,nx)
        END DO
     END IF
  END IF
  ! ... Precipitation flux
  ! ... KOD(18): 1 for spatially variable flux; 2 for uniform flux
  nprecip = 0
  IF (kod(18) > 0) THEN
     ! ... Set flags in IBC array for all the active precipitation nodes
     ! ...     Combination seepage and precipitation nodes already flagged 
     ! ...     by slice input of np
     DO j=1,ny
        DO i=1,nx
           ij = (j-1)*nx + i
           ikpt = (ktop(ij)-1)*nx + i       ! ... ik of top active cell
           IF(qprecip(ij) > 0._kdp) then
              ! ... Do not flag as precipitation flux 1)specified pressure cells
              ! ...      2) seepage face cells or 3) combination seepage and precipitation
              ! ...      cells
              if(ibc(ikpt,j) /= 11 .AND. ibc(ikpt,j) /= 30 .AND. ibc(ikpt,j) /= 50)  &
                   ibc(ikpt,j) = 20
              nprecip = nprecip + 1
           END IF
        END DO
     END DO
     ! ... Calculate specific enthalpy and density of precipitation for each
     ! ...      cell at top of region
     DO  ij=1,nxy
        CALL tempress(tflux(ij),patm,denflux(ij),ehflux(ij))
     END DO
     ! ... Check for potential ponding
     DO  ij=1,nxy
        CALL mtoijk(ij,i,j,k,nx,ny)
        icpt = (j-1)*nxz + (ktop(ij)-1)*nx + i
        visw(icpt) = 0.01_kdp     ! ... Preliminary set of viscosity to 1 P for check purposes
        IF(qprecip(ij) > zk(icpt)*denflux(ij)*grav/visw(icpt)) THEN
           WRITE(fustdout,9750) 'Precipitation flux too large, ',  &
                'ponding problem; i,j,k: ',i,j,ktop(ij),  &
                'Precipitation flux: ',qprecip(ij)
           WRITE(fuclog,9750) 'Precipitation flux too large, ',  &
                'ponding problem; i,j,k: ',i,j,ktop(ij),  &
                'Precipitation flux: ',qprecip(ij)
9750       FORMAT(tr5,2a,3i4,tr2,a,1pg14.6/)
           ierr(68) = .TRUE.
        END IF
     END DO
     IF (initial) THEN
        WRITE(fupdef,9005) 'Precipitation Flux', unitlabel(18)
     ELSE
        WRITE(fupdef,9005) 'New values for Precipitation Flux', unitlabel(18)
     END IF
     ! ... convert to user output units and write array
     IF (kod(18) == 2) THEN        ! ... uniform flux
        WRITE(fupdef,9025) qprecip(1)/unitfac(18)
     ELSE                          ! ... spatial distribution
        DO  j = 1,ny
           WRITE(fupdef, '(tr1)')
           WRITE(fupdef,9030) '--- Slice',j,' ---'
           ijxf = (j-1)*nx + 1
           ijxl = j*nx
           WRITE(fupdef,cfmt) (qprecip(ij)/unitfac(18),ij=ijxf,ijxl)
        END DO
     END IF
     IF (irad) THEN        ! ... top surface area for cylindrical coordinates
        WRITE(fupdef,'(/tr10,A/tr15,A)')  &
             'Top Surface Area per Cell for Cylindrical Coordinates (cm^2)',  &
             ' (Calculated from squared nodal cell boundaries)'
        WRITE(fupdef,'(1X)')
        WRITE(fupdef, cfmt) (pi*(rsq(i+1)-rsq(i)), i=1,nx)
     ELSE                  ! ... top surface area for Cartesian coordinates
        WRITE(fupdef, '(/tr10,A)')  &
             'Top Surface Area per Cell for Cartesian Coordinates (cm^2)'
        IF (ABS(kod(12)) == 2 .AND. ABS(kod(13)) == 2) THEN
           WRITE(fupdef,9025) dx(1)*dy(1)
        ELSE
           DO  j=1,ny
              WRITE(fupdef, '(tr1)')
              WRITE(fupdef,9030) '--- Slice',j,' ---'
              WRITE(fupdef, cfmt) (dx(i)*dy(j), i=1,nx)
           END DO
        END IF
     END IF
     ! ... Fluid input rate per node at top of active mesh
     WRITE(fupdef,'(/tr10,A)') 'Fluid Input Flow Rate per cell', '(cm^3/s)'
     IF (irad) THEN
        WRITE(fupdef,'(tr1)')
        WRITE(fupdef,cfmt) (qprecip(i)*pi*(rsq(i+1)-rsq(i)), i=1,nx)
     ELSE IF (kod(18) == 2 .AND. ABS(kod(12)) == 2 .AND. ABS(kod(13)) == 2) THEN
        WRITE(fupdef,9025) qprecip(1)*dx(1)*dy(1)
     ELSE               ! ... spatially variable
        DO  j=1,ny
           WRITE(fupdef,'(tr1)')
           WRITE(fupdef,9030) '--- Slice',j,' ---'
           WRITE(fupdef,cfmt) (qprecip((j-1)*nx+i)*dx(i)*dy(j), i=1,nx)
        END DO
     END IF
     ! ... Precipitation heat input rate per cell
     WRITE(fupdef,'(/tr10,A)') 'Precipitation Heat Input flow rate per cell (erg/s)'
     IF (irad) THEN
        WRITE(fupdef,'(tr1)')
        WRITE(fupdef,cfmt) (qprecip(i)*denflux(i)*ehflux(i)*  &
             pi*(rsq(i+1)-rsq(i)),i=1,nx)
     ELSE IF (kod(18) == 2 .AND. ABS(kod(12)) == 2 .AND. ABS(kod(13)) == 2) THEN
        WRITE(fupdef,9025) qprecip(1)*denflux(1)*ehflux(1)*dx(1)*dy(1)
     ELSE               ! ... spatial distribution
        DO  j = 1,ny
           WRITE(fupdef,'(tr1)')
           WRITE(fupdef,9030) '--- Slice',j,' ---'
           WRITE(fupdef,cfmt) (qprecip((j-1)*nx+i)*denflux((j-1)*nx+i)*  &
                ehflux((j-1)*nx+i)*dx(i)*dy(j),i=1,nx)
        END DO
     END IF
  END IF
9005 FORMAT (/'--- ',a,' ---  ',a)
9025 FORMAT (tr10,'Uniform Value -  ',1pg11.4)
9030 FORMAT (tr10,a,i4,a)
END SUBROUTINE init_wr_bc
