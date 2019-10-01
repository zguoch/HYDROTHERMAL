SUBROUTINE pdata(flag)
  ! ...  Purpose: To print data arrays as requested
  USE machine_constants, ONLY: kdp
  USE math_constants
  USE f_units
  USE bc, ONLY: ibc, seeping, cdtn, denflux, ehassoc, ehflux, qprecip, nhcond, nprecip, nseep,  &
       nspecpbc, unconfined
  USE parameters          !!$*** for flint
  USE control
  USE fdeq, ONLY: resmaxm, resmas, imres, jmres, kmres,  &
       resmaxe, reseng, ieres, jeres, keres,  &
       residm, reside
  USE mesh
!!$  USE parameters
  USE source
  USE variables
  USE units
  IMPLICIT NONE
  INTEGER, INTENT(IN OUT) :: flag
  !
  CHARACTER(LEN=80) :: cfmta
  CHARACTER(LEN=80), PARAMETER :: cfmta1='(i4,1P,1X,g8.2,i4,i3,i4,1X,a5,5(g10.2))',  &
       cfmta2='(i5,1P,2X,g9.3,3I4,2X,a7,5(2X,g10.3))'
  CHARACTER(LEN=1) :: rxlbl
  INTEGER :: i, ic, icp, icxf, icxl, ij, ik, ikp, j, k, kdum, klayers,  &
       nnodes, nr
  LOGICAL :: lprint
  REAL(KIND=kdp) ::  dum1, gveqne, gveqnead, gveqnecd, gveqnm,  &
       qcdtndum, qqh, qsppbc, timedisply, wmf, zsum
  REAL(KIND=kdp) :: qprecipdum, qhdum
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: lprint1
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: hvrtavg, pvrtavg, swvrtavg, tcvrtavg
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: dumarry, aprint1, aprint2
  INTEGER :: alloc_stat
  INTERFACE
     SUBROUTINE printar(ndim,array,lprnt,fu,cnv,jfmt,nnoppr)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: ndim
       REAL(kind=kdp), DIMENSION(:), INTENT(IN) :: array
       INTEGER, DIMENSION(:,:), INTENT(IN) :: lprnt
       INTEGER, INTENT(IN) :: fu
       REAL(KIND=kdp), INTENT(IN) :: cnv
       INTEGER, INTENT(IN) :: jfmt
       INTEGER, INTENT(IN) :: nnoppr
     END SUBROUTINE printar
  END INTERFACE     !!$*** for flint
  INTERFACE         !!$*** for flint
     SUBROUTINE wrnew(funit,dumarray,kodd,convfac)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: funit
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: dumarray
       INTEGER, INTENT(IN) :: kodd
       REAL(KIND=kdp), INTENT(IN OUT) :: convfac
     END SUBROUTINE wrnew
  END INTERFACE
  !     ------------------------------------------------------------------
  !       Ioptprta() - print option for arrays; 0 - skip; 1 - print
  !           (1)-pressure, (2)-enthalpy, (3)-temperature,
  !           (4)-water sat'n, (5)-resid mass, (6)-resid engy
  !           (7)-water density,  (8)-water viscosity,  (9)-water X velocity,
  !           (10)-water Y velocity, (11)-water Z velocity,
  !           (12)-water hydraulic potential
  !           (13)-steam density, (14)-steam viscosity, (15)-steam X velocity,
  !           (16)-steam Y velocity, (17)-steam Z velocity,
  !           (18)-steam hydraulic potential
  !           (19)-porosity, (20)-X-permeability, (21)-Y-permeability,
  !           (22)-Z-permeability, (23)-thermal conductivity,
  !           (24)-specific heat, (25)-rock density (26)-compressibility
  !           (27)-Cell Peclet number, (28)-Cell Nusselt number
  !
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.2 $//$Date: 2007/10/16 21:17:22 $'
  ALLOCATE(lprint1(nxz,ny), hvrtavg(nxy), pvrtavg(nxy), swvrtavg(nxy), tcvrtavg(nxy), &
           dumarry(nxyz), aprint1(nxyz), aprint2(nxyz), STAT=alloc_stat)
  IF (alloc_stat /= 0) then
     PRINT *, 'ARRAY ALLOCATION FAILED: pdata'
     ierr(188) = .TRUE.
     RETURN
  END IF
  !     ------------------------------------------------------------------
  !...
  ! ... Set value for display of time in seconds or years
  IF (iyr == 2) THEN
     timedisply = year(time)
  ELSE
     timedisply = time
  END IF
  IF(.NOT.irad) THEN
     rxlbl = 'X'
  ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
     nr = nx
     rxlbl = 'R'
  END IF
  ! ... assign formats for given printing width
  IF (ioptpr(2) == 0) THEN
     cfmta = cfmta1   ! ...     80 columns
  ELSE
     cfmta = cfmta2   ! ... 132 or 159 columns
  END IF
  ! ... Compute vertically averaged values
  IF (ioptpr(4) == 1) THEN
     DO  i=1,nx
        DO  j=1,ny
           zsum = 0._kdp
           ij = (j-1)*nx + i
           pvrtavg(ij) = 0.0_kdp
           hvrtavg(ij) = 0.0_kdp
           swvrtavg(ij) = 0.0_kdp
           tcvrtavg(ij) = 0.0_kdp
           DO  k=1,nz
              ik = (k-1)*nx + i
              IF (np(ik,j) /= 0) THEN
                 ic = (j-1)*nxz + (k-1)*nx + i
                 zsum = zsum + dz(k)
                 pvrtavg(ij) = pvrtavg(ij) + p(ic)*dz(k)
                 hvrtavg(ij) = hvrtavg(ij) + h(ic)*dz(k)
                 swvrtavg(ij) = swvrtavg(ij) + satnw(ic)*dz(k)
                 tcvrtavg(ij) = tcvrtavg(ij) + tc(ic)*dz(k)
              END IF
           END DO
           IF (zsum > 0._kdp) THEN
              pvrtavg(ij) = pvrtavg(ij)/zsum
              hvrtavg(ij) = hvrtavg(ij)/zsum
              swvrtavg(ij) = swvrtavg(ij)/zsum
              tcvrtavg(ij) = tcvrtavg(ij)/zsum
           END IF
        END DO
     END DO
  END IF
  lprint = .FALSE.
  ! ... Pressure
  IF(ABS(print_press) > 0._kdp) THEN
     IF(flag > 0) THEN
        lprint = .TRUE.
     ELSE
        CALL print_control(print_press,time,ltime,tchg,timprp,lprint)
     END IF
     IF (lprint .AND. ioptprta(1) == 1) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fup,2005) rxlbl//'-Direction Node Coordinates   '//TRIM(unitlabel(12))
2005          FORMAT(//tr30,A)
              CALL printar(1,xx,ibc,fup,unitfac(12),12,nx)
              WRITE(fup,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fup,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fup,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fup,unitfac(12),12,nr)
           END IF
           WRITE(fup,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fup,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fup,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fup,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        WRITE(fup,9010) 'Pressure Values', unitlabel(11)
        CALL printar(2,p,ibc,fup,unitfac(11),24,000)
        IF (ioptpr(4) == 1) THEN
           WRITE(fup,9010) 'Vertically Averaged Pressure', '(dyne/cm^2)'
           DO  j=1,ny
              icxf = (j-1)*nx + 1
              icxl = icxf + nx - 1
              WRITE(fup,9005) (pvrtavg(ij), ij=icxf, icxl)
           END DO
        END IF
        ! ... Water-table elevation
        IF(unconfined) THEN
           CALL wt_elev
           lprint1 = -1
           DO  j = 1,ny
              DO  k = nz,1,-1
                 DO  i = 1,nx
                    ik = (k-1)*nx + i
                    ic = (j-1)*nxz + (k-1)*nx + i
                    IF (ABS(fs_elev(ic)) > 0._kdp) THEN
                       lprint1(ik,j) = 1
                    END IF
                 END DO
              END DO
           END DO
           IF (iyr == 2) THEN
              WRITE(fuwtelev,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
           ELSE
              WRITE(fuwtelev,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
           END IF
           WRITE(fuwtelev,9010) 'Water-Table Elevation Values', unitlabel(14)
           CALL printar(2,fs_elev,lprint1,fuwtelev,unitfac(14),24,000)
           IF(ltime > 0) THEN
              WRITE(fustdout,9155) ltime, 'L2-norm of change in water-table elevation:',  &
                   dfs_l2norm/unitfac(14), unitlabel(14)
              WRITE(fuclog,9155) ltime,  'L2-norm of change in water-table elevation:',  &
                   dfs_l2norm/unitfac(14), unitlabel(14)
9155          FORMAT(i4,tr4,a,1pg15.4,a)
           END IF
        END IF
     END IF
  END IF
  ! ... Enthalpy
  IF(ABS(print_enth) > 0.) THEN
     IF(flag > 0) THEN
        lprint = .TRUE.
     ELSE
        CALL print_control(print_enth,time,ltime,tchg,timprh,lprint)
     END IF
     IF (lprint .AND. ioptprta(2) == 1) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fuen,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fuen,unitfac(12),12,nx)
              WRITE(fuen,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fuen,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fuen,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fuen,unitfac(12),12,nr)
           END IF
           WRITE(fuen,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fuen,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fuen,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fuen,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        WRITE(fuen,9010) 'Enthalpy Values', unitlabel(10)
        CALL printar(2,h,ibc,fuen,unitfac(10),24,000)
        IF (ioptpr(4) == 1) THEN
           ! ... Print vertically averaged enthalpy
           WRITE(fuen,9010) 'Vertically Averaged Enthalpy', '(erg/g)'
           DO  j=1,ny
              icxf = (j-1)*nx + 1
              icxl = icxf + nx - 1
              WRITE(fuen,9005) (hvrtavg(ij), ij=icxf, icxl)
           END DO
        END IF
     END IF
  END IF
  ! ... Temperature
  IF(ABS(print_temp) > 0.) THEN
     CALL print_control(print_temp,time,ltime,tchg,timprt,lprint)
     IF(flag > 0) lprint = .TRUE.
     IF (lprint .AND. ioptprta(3) == 1) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fut,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fut,unitfac(12),12,nx)
              WRITE(fut,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fut,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fut,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fut,unitfac(12),12,nr)
           END IF
           WRITE(fut,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fut,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fut,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fut,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        WRITE(fut,9010) 'Temperature Values', unitlabel(15)
        CALL printar(2,tc,ibc,fut,unitfac(15),11,000)
        IF (ioptpr(4) == 1) THEN
           ! ... Print vertically averaged temperature
           WRITE(fut,9010) 'Vertically Averaged Temperature', '(Deg.C)'
           DO  j=1,ny
              icxf = (j-1)*nx + 1
              icxl = icxf + nx - 1
              WRITE(fut,9005) (tcvrtavg(ij), ij=icxf, icxl)
           END DO
        END IF
     END IF
  END IF
  ! ... Saturation
  IF(ABS(print_satn) > 0.) THEN
     CALL print_control(print_satn,time,ltime,tchg,timprsat,lprint)
     IF(flag > 0) lprint = .TRUE.
     IF (lprint .AND. ioptprta(4) > 0) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fusat,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fusat,unitfac(12),12,nx)
              WRITE(fusat,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fusat,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fusat,2005) rxlbl//'-Direction Node Coordinate   '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fusat,unitfac(12),12,nr)
           END IF
           WRITE(fusat,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fusat,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fusat,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fusat,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        IF (ioptprta(4) == 1 .OR. ioptprta(4) == 3) THEN        ! ... water saturation (volumetric)
           WRITE(fusat,9010) 'Water Saturation Values', '(volumetric fraction)'
           CALL printar(2,satnw,ibc,fusat,1._kdp,14,000)
           IF (ioptpr(4) == 1) THEN          ! ... Print vertically averaged water saturation
              WRITE(fusat,9010) 'Vertically Averaged Water Saturation','(volumetric fraction)'
              DO  j=1,ny
                 icxf = (j-1)*nx + 1
                 icxl = icxf + nx - 1
                 WRITE(fusat,9005) (swvrtavg(ij), ij=icxf, icxl)
              END DO
           END IF
        END IF
        IF (ioptprta(4) == 2 .OR. ioptprta(4) == 3) THEN     ! ... Mass fraction water
           DO  ic=1,nxyz
              dumarry(ic) = satnw(ic)*denw(ic)/  &
                   (satnw(ic)*denw(ic)+(1.0_KDP-satnw(ic))*dens(ic))
           END DO
           WRITE(fusat,9010) 'Mass Fraction Water', '(-)'
           CALL printar(2,dumarry,ibc,fusat,1._kdp,24,000)
        END IF
     END IF
  END IF
  ! ... Residuals
  IF(ABS(print_resid) > 0.) THEN
     CALL print_control(print_resid,time,ltime,tchg,timprresid,lprint)
     IF(flag > 0) lprint = .TRUE.
     IF (lprint) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fures,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fures,unitfac(12),12,nx)
              WRITE(fures,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fures,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fures,2005) rxlbl//'-Direction Node Coordinate '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fures,unitfac(12),12,nr)
           END IF
           WRITE(fures,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fures,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fures,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fures,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        IF (ioptprta(5) == 1) THEN       ! ... Residual mass balance
           IF (ltime == 0) THEN
              WRITE(fures, '(10X,A)')  &
                   'Residuals are not calculated for the initial values'
           ELSE
              WRITE(fures,9010) 'Residual Mass Balance', '(g/s-cm^3)'
              WRITE(fures,9020) resmaxm, imres, jmres, kmres, resmas
              CALL printar(2,residm,ibc,fures,1._kdp,24,000)
           END IF
        END IF
        IF (ioptprta(6) == 1) THEN        ! ... Residual energy balance
           IF (ltime == 0) THEN
              WRITE(fures, '(10X,A)')  &
                   'Residuals are not calculated for the initial values'
           ELSE
              WRITE(fures,9010) 'Residual Energy Balance', '(erg/s-cm^3)'
              WRITE(fures,9020) resmaxe, ieres, jeres, keres, reseng
              CALL printar(2,reside,ibc,fures,1._kdp,24,000)
           END IF
        END IF
     END IF
  END IF
  ! ... Density of water and steam
  IF(ABS(print_dens) > 0.) THEN
     CALL print_control(print_dens,time,ltime,tchg,timprd,lprint)
     IF(flag > 0) lprint = .TRUE.
     IF (lprint) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fuden,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fuden,unitfac(12),12,nx)
              WRITE(fuden,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fuden,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fuden,2005) rxlbl//'-Direction Node Coordinate '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fuden,unitfac(12),12,nr)
           END IF
           WRITE(fuden,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fuden,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fuden,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fuden,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        IF (ioptprta(7) == 1) THEN
           WRITE(fuden,9010) 'Water Density Values', '(g/cm^3)'
           CALL printar(2,denw,ibc,fuden,1._kdp,14,000)
        END IF
        IF (ioptprta(13) == 1) THEN
           WRITE(fuden,9010) 'Steam Density Values', '(g/cm^3)'
           CALL printar(2,dens,ibc,fuden,1._kdp,14,000)
        END IF
     END IF
  END IF
  ! ... Viscosity of water and steam
  IF(ABS(print_vis) > 0.) THEN
     CALL print_control(print_vis,time,ltime,tchg,timprv,lprint)
     IF(flag > 0) lprint = .TRUE.
     IF (lprint) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fuvis,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fuvis,unitfac(12),12,nx)
              WRITE(fuvis,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fuvis,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fuvis,2005) rxlbl//'-Direction Node Coordinate '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fuvis,unitfac(12),12,nr)
           END IF
           WRITE(fuvis,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fuvis,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fuvis,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fuvis,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        IF (ioptprta(8) == 1) THEN
           WRITE(fuvis,9010) 'Water Viscosity Values', '(poise);(g/cm-s)'
           CALL printar(2,visw,ibc,fuvis,1._kdp,24,000)
        END IF
        IF (ioptprta(14) == 1) THEN
           WRITE(fuvis,9010) 'Steam Viscosity Values', '(poise);(g/cm-s)'
           CALL printar(2,viss,ibc,fuvis,1._kdp,24,000)
        END IF
     END IF
  END IF
  ! ... Velocity of water and steam
  IF(ABS(print_vel) > 0.) THEN
     CALL print_control(print_vel,time,ltime,tchg,timprvel,lprint)
     IF(flag > 0) lprint = .TRUE.
     !..            if(nx.eq.1) ioptprta(9) = 0
     IF (lprint) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fuvel,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fuvel,unitfac(12),12,nx)
              WRITE(fuvel,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fuvel,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fuvel,2005) rxlbl//'-Direction Node Coordinate '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fuvel,unitfac(12),12,nr)
           END IF
           WRITE(fuvel,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fuvel,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fuvel,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fuvel,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        IF ((ioptprta(9) == 1 .OR. ioptprta(9) == 3)) THEN     ! ... X-velocity of water
           WRITE(fuvel,9010) 'X Water Interstitial Velocity', '(cm/s)'
           WRITE(fuvel,9035) xwvelmax
           CALL printar(2,xwvel,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(10) == 1 .OR. ioptprta(10) == 3)) THEN        ! ... Y-velocity of water
           WRITE(fuvel,9010) 'Y Water Interstitial Velocity', '(cm/s)'
           WRITE(fuvel,9035) ywvelmax
           CALL printar(2,ywvel,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(11) == 1 .OR. ioptprta(11) == 3)) THEN     ! ... Z-velocity of water
           WRITE(fuvel,9010) 'Z Water Interstitial Velocity (-=down, +=up)', '(cm/s)'
           WRITE(fuvel,9035) zwvelmax
           CALL printar(2,zwvel,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(9) == 2 .OR. ioptprta(9) == 3)) THEN     ! ... X mass flux of water
           WRITE(fuvel,9010) 'X Water Mass Flux', '(g/s-cm^2)'
           WRITE(fuvel,9035) xwmflxmax
           CALL printar(2,xwmflx,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(10) == 2 .OR. ioptprta(10) == 3)) THEN     ! ... Y mass flux of water
           WRITE(fuvel,9010) 'Y Water Mass Flux', '(g/s-cm^2)'
           WRITE(fuvel,9035) ywmflxmax
           CALL printar(2,ywmflx,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(11) == 2 .OR. ioptprta(11) == 3)) THEN     ! ... Z mass flux of water
           WRITE(fuvel,9010) 'Z Water Mass Flux (-=down, +=up)', '(g/s-cm^2)'
           WRITE(fuvel,9035) zwmflxmax
           CALL printar(2,zwmflx,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(9) == 2 .OR. ioptprta(9) == 3)) THEN       ! ... X Darcy flux of water 
           WRITE(fuvel,9010) 'X Water Darcy Flux at cell face between x(i-1) and x(i)'//  &
                ' (cm^3/s-cm^2)'
           CALL printar(2,xwflux,lprntxif,fuvel,1._kdp,24,200)
        END IF
        IF ((ioptprta(10) == 2 .OR. ioptprta(10) == 3)) THEN        ! ... Y Darcy flux of water
           WRITE(fuvel,9010) 'Y Water Darcy Flux at cell face between y(j-1) and y(j)'//  &
                ' (cm^3/s-cm^2)'
           CALL printar(2,ywflux,lprntyif,fuvel,1._kdp,24,020)
        END IF
        IF ((ioptprta(11) == 2 .OR. ioptprta(11) == 3)) THEN     ! ... Z Darcy flux of water
           WRITE(fuvel,9010) 'Z Water Darcy Flux at cell face between z(k-1) and z(k)'//  &
                ' (-=down, +=up)', '(cm^3/s-cm^2)'
           CALL printar(2,zwflux,lprntzif,fuvel,1._kdp,24,002)
        END IF
        IF ((ioptprta(15) == 1 .OR. ioptprta(15) == 3)) THEN     ! ... X-velocity of steam
           WRITE(fuvel,9010) 'X Steam Interstitial Velocity', '(cm/s)'
           WRITE(fuvel,9035) xsvelmax
           CALL printar(2,xsvel,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(16) == 1 .OR. ioptprta(16) == 3)) THEN     ! ... Y-velocity of steam
           WRITE(fuvel,9010) 'Y Steam Interstitial Velocity', '(cm/s)'
           WRITE(fuvel,9035) ysvelmax
           CALL printar(2,ysvel,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(17) == 1 .OR. ioptprta(17) == 3)) THEN     ! ... Z-velocity of steam
           WRITE(fuvel,9010) 'Z Steam Interstitial Velocity (-=down, +=up)', '(cm/s)'
           WRITE(fuvel,9035) zsvelmax
           CALL printar(2,zsvel,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(15) == 2 .OR. ioptprta(15) == 3)) THEN     ! ... X mass flux of steam
           WRITE(fuvel,9010) 'X Steam Mass Flux', '(g/s-cm^2)'
           WRITE(fuvel,9035) xsmflxmax
           CALL printar(2,xsmflx,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(16) == 2 .OR. ioptprta(16) == 3)) THEN     ! ... Y mass flux of steam
           WRITE(fuvel,9010) 'Y Steam Mass Flux', '(g/s-cm^2)'
           WRITE(fuvel,9035) ysmflxmax
           CALL printar(2,ysmflx,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(17) == 2 .OR. ioptprta(17) == 3)) THEN     ! ... Z mass flux of steam
           WRITE(fuvel,9010) 'Z Steam Mass Flux (-=down, +=up)', '(g/s-cm^2)'
           WRITE(fuvel,9035) zsmflxmax
           CALL printar(2,zsmflx,ibc,fuvel,1._kdp,24,000)
        END IF
        IF ((ioptprta(15) == 2 .OR. ioptprta(15) == 3)) THEN       ! ... X Darcy flux of steam
           WRITE(fuvel,9010) 'X Steam Darcy Flux at cell face between x(i-1) and x(i)'//  &
                ' (cm^3/s-cm^2)'
           CALL printar(2,xsflux,lprntxif,fuvel,1._kdp,24,200)
        END IF
        IF ((ioptprta(16) == 2 .OR. ioptprta(16) == 3)) THEN        ! ... Y Darcy flux of steam
           WRITE(fuvel,9010) 'Y Steam Darcy Flux at cell face between y(j-1) and y(j)'//  &
                ' (cm^3/s-cm^2)'
           CALL printar(2,ysflux,lprntyif,fuvel,1._kdp,24,020)
        END IF
        IF ((ioptprta(17) == 2 .OR. ioptprta(17) == 3)) THEN     ! ... Z Darcy flux of steam
           WRITE(fuvel,9010) 'Z Steam Darcy Flux at cell face between z(k-1) and z(k)'//  &
                ' (-=down, +=up)', '(cm^3/s-cm^2)'
           CALL printar(2,zsflux,lprntzif,fuvel,1._kdp,24,002)
        END IF
     END IF
  END IF
  ! ... Hydraulic potential of water and steam
  IF(ABS(print_pot) > 0._kdp) THEN
     CALL print_control(print_pot,time,ltime,tchg,timprpot,lprint)
     IF(flag > 0) lprint = .TRUE.
     IF(lprint) CALL hyd_potential     ! ... calculate the hydraulic potential of steam and water
     IF (lprint) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fupot,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fupot,unitfac(12),12,nx)
              WRITE(fupot,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fupot,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fupot,2005) rxlbl//'-Direction Node Coordinate '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fupot,unitfac(12),12,nr)
           END IF
           WRITE(fupot,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fupot,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fupot,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fupot,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        IF (ioptprta(12) == 1 .OR. ioptprta(12) == 3) THEN
           WRITE(fupot,9010) 'Water Hydraulic Potential Values',unitlabel(11)
           CALL printar(2,wpot,ibc,fupot,unitfac(11),24,000)
        END IF
        IF (ioptprta(12) == 2 .OR. ioptprta(12) == 3) THEN
           WRITE(fupot,9010) 'Water Potentiometric Head Values',unitlabel(14)
           CALL printar(2,hwpot,ibc,fupot,unitfac(14),24,000)
        END IF
        IF (ioptprta(18) == 1 .OR. ioptprta(18) == 3) THEN
           WRITE(fupot,9010) 'Steam Hydraulic Potential Values',unitlabel(11)
           CALL printar(2,spot,ibc,fupot,unitfac(11),24,000)
        END IF
        IF (ioptprta(18) == 2 .OR. ioptprta(18) == 3) THEN
           WRITE(fupot,9010) 'Steam Potentiometric Head Values',unitlabel(14)
           CALL printar(2,hspot,ibc,fupot,unitfac(14),24,000)
        END IF
     END IF
  END IF
  ! ... Porosity
  IF(ABS(print_poros) > 0.) THEN
     CALL print_control(print_poros,time,ltime,tchg,timprpor,lprint)
     IF(flag > 0) lprint = .TRUE.
     IF (lprint) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fupor,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fupor,unitfac(12),12,nx)
              WRITE(fupor,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fupor,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fupor,2005) rxlbl//'-Direction Node Coordinate '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fupor,unitfac(12),12,nr)
           END IF
           WRITE(fupor,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fupor,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fupor,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fupor,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        IF (ioptprta(19) == 1) THEN
           WRITE(fupor,9010) 'Media Porosity', unitlabel(1)
           IF(kod(1) == 1) THEN
              CALL printar(2,phi,ibc,fupor,unitfac(1),14,000)
           ELSE
              CALL wrnew(fupor, phi, kod(1), unitfac(1))
           END IF
        END IF
     END IF
  END IF
  ! ... Permeability
  IF(ABS(print_perm) > 0.) THEN
     CALL print_control(print_perm,time,ltime,tchg,timprperm,lprint)
     IF(flag > 0) lprint = .TRUE.
     IF (lprint) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fuperm,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fuperm,unitfac(12),12,nx)
              WRITE(fuperm,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fuperm,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fuperm,2005) rxlbl//'-Direction Node Coordinate '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fuperm,unitfac(12),12,nr)
           END IF
           WRITE(fuperm,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fuperm,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fuperm,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fuperm,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        IF (ioptprta(20) >= 1) THEN     ! ... x-permeability
           WRITE(fuperm,9010) 'X - Permeability', unitlabel(2)
           IF(kod(2) == 1) THEN
              CALL printar(2,xk,ibc,fuperm,unitfac(2),24,000)
           ELSE
              CALL wrnew(fuperm, xk, kod(2), unitfac(2))
           END IF
        END IF
        IF (ioptprta(21) >= 1) THEN     ! ... y-permeability
           WRITE(fuperm,9010) 'Y - Permeability', unitlabel(3)
           IF(kod(3) == 1) THEN
              CALL printar(2,yk,ibc,fuperm,unitfac(3),24,000)
           ELSE
              CALL wrnew(fuperm, yk, kod(3), unitfac(3))
           END IF
        END IF
        IF (ioptprta(22) >= 1) THEN     ! ... z-permeability
           WRITE(fuperm,9010) 'Z - Permeability', unitlabel(4)
           IF(kod(4) == 1) THEN
              CALL printar(2,zk,ibc,fuperm,unitfac(4),24,000)
           ELSE
              CALL wrnew(fuperm, zk, kod(4), unitfac(4))
           END IF
        END IF
        ! ... Relative permeability 
        IF (ioptprta(20) > 0 .OR. ioptprta(21) > 0 .OR. ioptprta(22) > 0) THEN
           WRITE(fuperm,9010) 'Water Relative Permeability'
           CALL printar(2,rw,ibc,fuperm,1._kdp,24,000)
           WRITE(fuperm,9010) 'Steam Relative Permeability'
           CALL printar(2,rs,ibc,fuperm,1._kdp,24,000)
        END IF
     END IF
  END IF
  ! ... Porous media properties
  IF(ABS(print_pmprop) > 0.) THEN
     CALL print_control(print_pmprop,time,ltime,tchg,timprprop,lprint)
     IF(flag > 0) lprint = .TRUE.
     IF (lprint) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(futhp,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,futhp,unitfac(12),12,nx)
              WRITE(futhp,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,futhp,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(futhp,2005) rxlbl//'-Direction Node Coordinate '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,futhp,unitfac(12),12,nr)
           END IF
           WRITE(futhp,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,futhp,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(futhp,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(futhp,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        IF (ioptprta(23) == 1) THEN     ! ...      Thermal Conductivity
           WRITE(futhp,9010) 'Porous Media Thermal Conductivity', unitlabel(5)
           IF(kod(5) == 1) THEN
              CALL printar(2,xkc,ibc,futhp,unitfac(5),24,000)
           ELSE
              CALL wrnew(futhp, xkc, kod(5), unitfac(5))
           END IF
        END IF
        IF (ioptprta(24) == 1) THEN     ! ...  Specific Heat of Rock
           WRITE(futhp,9010) 'Specific Heat of Rock', unitlabel(6)
           IF(kod(6) == 1) THEN
              CALL printar(2,phfwt,ibc,futhp,unitfac(6),24,000)
           ELSE
              CALL wrnew(futhp, phfwt, kod(6), unitfac(6))
           END IF
        END IF
        IF (ioptprta(25) == 1) THEN     ! ...  Rock Density
           WRITE(futhp,9010) 'Rock Density', unitlabel(7)
           IF(kod(7) == 1) THEN
              CALL printar(2,df,ibc,futhp,unitfac(7),24,000)
           ELSE
              CALL wrnew(futhp, df, kod(7), unitfac(7))
           END IF
        END IF
        IF (ioptprta(26) == 1) THEN     ! ...  Matrix Compressibility
           WRITE(futhp,9010) 'Porous Matrix Compressibility', unitlabel(8)
           IF(kod(8) == 1) THEN
              CALL printar(2,beta,ibc,futhp,unitfac(8),24,000)
           ELSE
              CALL wrnew(futhp, beta, kod(8), unitfac(8))
           END IF
        END IF
     END IF
  END IF
  ! ...  Time series output
  IF(ABS(print_temporal) > 0.) THEN
     CALL print_control(print_temporal,time,ltime,tchg,timprtemp,lprint)
     IF(flag > 0) lprint = .TRUE.      ! ... force output if error exit
     IF (lprint) THEN
        IF (ioptpr(6) > 0) THEN        ! ... For selected nodes
           DO  kdum = 1,ioptpr(6)
              i = ijktp(kdum,1)
              j = ijktp(kdum,2)
              k = ijktp(kdum,3)
              ic = (k-1)*nx + i + (j-1)*nxz
              ! ... Compute water mass fraction
              wmf = satnw(ic)*denw(ic)/(satnw(ic)*denw(ic)+(1.0_kdp-satnw(ic))*dens(ic))
              IF (ioptpr(2) >= 1) THEN
                 ! ... for wide (132 column) output
                 WRITE(fuplttimser,9030) ltime, timedisply, i, j, k,  &
                      p(ic)/unitfac(11), h(ic)/unitfac(10),  &
                      tc(ic)/unitfac(15), satnw(ic), wmf
9030             FORMAT(i5,1pe11.4,3I3,2(tr2,1pg11.4),0pf10.2,2F13.6)
              ELSE
                 ! ... for narrow (80 column) output
                 WRITE(fuplttimser,9025) ltime, timedisply, i, j, k,  &
                      p(ic)/unitfac(11), h(ic)/unitfac(10), tc(ic)/unitfac(15), satnw(ic)
9025             FORMAT(i5,1pe11.4,3I3,2(tr2,1pg11.4),0pf10.2,0pf13.6)
              END IF
           END DO
        ELSEIF (ioptpr(6) == -1) THEN        ! ...  For all nodes
           DO  j = 1, ny
              DO  k = 1, nz
                 DO  i = 1, nx
                    ic = (k-1)*nx + i + (j-1)*nxz
                    ! ... Compute water mass fraction
                    wmf = satnw(ic)*denw(ic)/(satnw(ic)*denw(ic)+(1.0_kdp-satnw(ic))*dens(ic))
                    IF (ioptpr(2) >= 1) THEN
                       ! ... For wide (132 column) output
                       WRITE(fuplttimser,9030) ltime, timedisply, i, j, k,  &
                            p(ic)/unitfac(11), h(ic)/unitfac(10),  &
                            tc(ic)/unitfac(15), satnw(ic), wmf
                    ELSE
                       ! ... For narrow (80 column) output
                       WRITE(fuplttimser,9025) ltime, timedisply, i, j, k,  &
                            p(ic)/unitfac(11), h(ic)/unitfac(10), tc(ic)/unitfac(15), satnw(ic)
                    END IF
                 END DO
              END DO
           END DO
        END IF
     END IF
  END IF
  ! ... Flux data for source/sink & constant P & H nodes and seepage nodes
  IF(ABS(print_source) > 0.) THEN
     CALL print_control(print_source,time,ltime,tchg,timprsrc,lprint)
     IF(flag > 0) lprint = .TRUE.
  END IF
  !*****Need to revise
  !*****Logical difficulty with sources and fluxes stored in same data structure
  !*****     must specify both types of printout to get either one
  IF (lnri >= 1 .AND. (lprint .AND. ioptpr(5) == 1)) THEN
     WRITE(fusrc,'(1x)')
     ! ... Print source/sink data for individual nodes and wells
     DO  j = 1,ny
        DO  i = 1,nx
           ij = (j-1)*nx + i
           IF(ABS(qt(ij)) > 0._kdp) THEN
              ! ...   Count the number of open intervals for this well
              klayers = 0
              DO  k = 1,nz
                 ic = (k-1)*nx + i + (j-1)*nxz
                 IF (iq(ic) == 1) klayers = klayers + 1
              END DO
           END IF
           IF (qt(ij) > 0.0_kdp) THEN              ! ... Injection well
              ! ... if more than one open interval in well, write values for total injection
              IF (klayers > 1) WRITE(fusrc,cfmta) ltime, timedisply, i,  &
                   j, ktop(ij), 'Injec Well', qt(ij), qhi(ij)*qt(ij), qhi(ij)*qt(ij),  &
                   0.0_kdp, qhi(ij)
           ELSEIF (qt(ij) < 0.0_kdp) THEN          ! ... Production well
              ! ... sum the layer flow rates
              qqh = 0.0_kdp
              DO  k = 1, nz
                 ic = (k-1)*nx + i + (j-1)*nxz
                 qqh = qqh + qh(ic)
              END DO
              ! ... if more than one open interval in well, write values for total production
              IF (klayers > 1) WRITE(fusrc,cfmta) ltime, timedisply, i,  &
                   j, ktop(ij), 'Prod Well', qt(ij), qqh, qqh, 0.0_kdp, qqh/qt(ij)
           END IF
           ! ... Individual source/sink nodes
           ! ...    for the bottom row, separate and print the conductive
           ! ...    heat flow through the boundary
           k = 1     ! ... Bottom row
           ic = (k-1)*nx + i + (j-1)*nxz
           IF(ABS(q(ic)) > 0.0_kdp) THEN
              ! ... Calculate conductive heat component
              IF (irad) THEN
                 qcdtndum = cdtn(ij)*pi*(rsq(i+1)-rsq(i))
              ELSE
                 qcdtndum = cdtn(ij)*dx(i)*dy(j)
              END IF
              IF(q(ic) > 0.0_kdp) WRITE(fusrc,cfmta) ltime, timedisply,  &
                   i, j, k, 'Source ', q(ic), qh(ic), qh(ic) - qcdtndum,  &
                   qcdtndum, (qh(ic)-qcdtndum)/q(ic)
              IF(q(ic) < 0.0_kdp) WRITE(fusrc,cfmta) ltime, timedisply,  &
                   i, j, k, 'Sink   ', q(ic), qh(ic), qh(ic)+qcdtndum,  &
                   qcdtndum, (qh(ic)+qcdtndum)/q(ic)
           END IF
           DO k=nz-1,2,-1
              ! ... All rows except the top and bottom rows
              ic = (k-1)*nx + i + (j-1)*nxz
              IF (q(ic) > 0.0_kdp) WRITE(fusrc,cfmta) ltime, timedisply,  &
                   i, j, k, 'Source ', q(ic), qh(ic), qh(ic), 0.0_kdp,  &
                   qh(ic)/q(ic)
              IF (q(ic) < 0.0_kdp) WRITE(fusrc,cfmta) ltime, timedisply,  &
                   i, j, k, 'Sink   ', q(ic), qh(ic), qh(ic), 0.0_kdp,  &
                   qh(ic)/q(ic)
           END DO
           k = nz     ! ...  Top row, treat like internal point source
           ic = (k-1)*nx + i + (j-1)*nxz
           IF(ABS(q(ic)) > 0.0_kdp) THEN
!!$              !c.....calculate precipitation flow component and advective heat flow
!!$              !              IF (Irad) THEN
!!$              !                qprecipdum = Denflux(ij)*Qprecip(ij)*
!!$              !     &               pi*(Rsq(i+1)-Rsq(i))
!!$              !              ELSE
!!$              !                qprecipdum = Denflux(ij)*Qprecip(ij)*Dx(i)*Dy(j)
!!$              !              ENDIF
!!$              !                qhdum = Ehflux(ij)*qprecipdum
              IF(ABS(qprecip(ij)) == 0._kdp) THEN
                 IF (q(ic) > 0.0_kdp) WRITE(fusrc,cfmta) ltime, timedisply,  &
                      i, j, k, 'Source ', q(ic), qh(ic), qh(ic), 0.0_kdp,  &
                      qh(ic)/q(ic)
                 IF (q(ic) < 0.0_kdp) WRITE(fusrc,cfmta) ltime, timedisply,  &
                      i, j, k, 'Sink   ', q(ic), qh(ic), qh(ic), 0.0_kdp,  &
                      qh(ic)/q(ic)
                 !              IF (Q(ic).GT.0.0_KDP) WRITE(fusrc, cfmta) Ltime, timedisply,
                 !     &                                   i, j, k, 'Source ', Q(ic),
                 !     &                                   Q(ic) - qprecipdum,
                 !     &                                   qprecipdum, (Q(ic)-qprecipdum)
                 !     &                                   /Q(ic)
                 ! ... Evapotranspiration (not implemented yet)
              END IF
           END IF
        END DO
     END DO
  END IF
  IF(ABS(print_bcflow) > 0.) THEN
     CALL print_control(print_bcflow,time,ltime,tchg,timprbcf,lprint)
     IF(flag > 0) lprint = .TRUE.
     IF(lprint) THEN
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fubcf2,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fubcf2,unitfac(12),12,nx)
              WRITE(fubcf2,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fubcf2,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fubcf2,2005) rxlbl//'-Direction Node Coordinate '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fubcf2,unitfac(12),12,nr)
           END IF
           WRITE(fubcf2,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fubcf2,unitfac(14),12,nz)
        END IF
        WRITE(fubcf2,2001) '*** Output at End of Time Step No. ', ltime,' ***'
2001    FORMAT(/tr30,a,i5,a)
        IF (iyr == 2) THEN
           WRITE(fubcf2,2002) 'Time '//dots,year(time),'(yr)'
        ELSE
           WRITE(fubcf2,2002) 'Time '//dots,time,'(s)'
        END IF
2002    FORMAT(/tr5,a60,1pg12.3,tr1,a)
        IF(COUNT(ibc == 11) > 0) THEN        ! ... Constant value nodes
           aprint1 = 0._kdp
           aprint2 = 0._kdp
           lprint1 = -1
           DO  j = 1,ny
              DO  k = nz,1,-1
                 DO  i = 1,nx
                    ik = (k-1)*nx + i
                    ic = (j-1)*nxz + (k-1)*nx + i
                    ! ... select  constant P and H nodes
                    IF (ibc(ik,j) == 11) THEN
                       CALL goveqn2(gveqnm,gveqne,gveqnead,gveqnecd,nnodes,i,j,k)
                       gveqnm = -gveqnm
                       gveqne = -gveqne
                       gveqnead = -gveqnead
                       gveqnecd = -gveqnecd
                       aprint1(ic) = gveqnm
                       aprint2(ic) = gveqne
                       lprint1(ik,j) = 1
                       IF (gveqnm /= 0._kdp) THEN
                          !                  if the number of adjacent active nodes is >1 then
                          !                    write that number rather than the enthalpy
                          IF (nnodes >= 2) THEN
                             dum1 = DBLE(nnodes)
                          ELSE
                             dum1 = gveqnead/gveqnm
                          END IF
                          WRITE(fubcf, cfmta) ltime, timedisply, i, j, k,  &
                               'ConsP&H', gveqnm, gveqne, gveqnead, gveqnecd, dum1
                       ELSE
                          WRITE(fubcf, cfmta) ltime, timedisply, i, j, k,  &
                               'ConsP&H', gveqnm, gveqne, gveqnead, gveqnecd
                       END IF
                    END IF
                 END DO
              END DO
           END DO
           WRITE(fubcf2,9110) 'Specified Pressure and Enthalpy,Temperature ',  &
                'B.C. Flow Rates (average over time step)','(positive is into the region)'
9110       FORMAT(//tr25,2a/tr50,a)
           WRITE(fubcf2,2042) 'Fluid   (g/s)'
2042       FORMAT(tr20,10A)
           CALL printar(2,aprint1,lprint1,fubcf2,1._kdp,24,000)
           WRITE(fubcf2,2042) 'Heat   (erg/s)'
           CALL printar(2,aprint2,lprint1,fubcf2,1._kdp,24,000)
        END IF
        IF(nspecpbc > 0) THEN     ! ... Specified pressure b.c. nodes with associated 
           !                             enthalpy, temperature
           aprint1 = 0._kdp
           aprint2 = 0._kdp
           lprint1 = -1
           DO  j = 1,ny
              DO  k = nz,1,-1
                 DO  i = 1,nx
                    ik = (k-1)*nx + i
                    ic = (j-1)*nxz + (k-1)*nx + i
                    ! ... Select specifed pressure nodes
                    IF (ibc(ik,j) == 10) THEN
                       CALL gov_spp(gveqnm,i,j,k)
                       ic = (k-1)*nx + i + (j-1)*nxz
                       ! ... Add in the capacitance and source/sink terms
                       gveqnm = gveqnm - (xm(ic)-xmoldt(ic))/delt + q(ic)
                       qsppbc = -gveqnm      ! ... Change sign so inflow is positive
                       IF(qsppbc > 0._kdp) THEN     ! ... inflow
                          gveqne = qsppbc*ehassoc(ic)
                       ELSE                         ! ... outflow
                          gveqne = qsppbc*h(ic)
                       END IF
                       aprint1(ic) = qsppbc
                       aprint2(ic) = gveqne
                       lprint1(ik,j) = 1
                       IF (ioptpr(2) == 0) THEN  ! ... 80 columns
                          WRITE(fubcf,9950) ltime, timedisply, i, j, k,  &
                              'SpPasoH', qsppbc, gveqne
                       ELSE                    ! ... 132 or 159 columns
                          WRITE(fubcf,9955) ltime, timedisply, i, j, k,  &
                              'SpPasoH', qsppbc, gveqne
                       END IF
                    END IF
                 END DO
              END DO
           END DO
           WRITE(fubcf2,2043)  &
                'Specified Pressure B.C. Flow Rates (at end of time step)',  &
                '(positive is into the region)'
           WRITE(fubcf2,2042) 'Fluid Mass   (g/s)'
           CALL printar(2,aprint1,lprint1,fubcf2,1._kdp,24,000)
           WRITE(fubcf2,2042) 'Associated Advective Heat   (erg/s)'
           CALL printar(2,aprint2,lprint1,fubcf2,1._kdp,24,000)
        END IF
        IF(nseep > 0) THEN        ! ... Seepage surface nodes
           aprint1 = 0._kdp
           aprint2 = 0._kdp
           lprint1 = -1
           DO  j = 1,ny
              DO  k = nz,1,-1
                 DO  i = 1,nx
                    ik = (k-1)*nx + i
                    ic = (j-1)*nxz + (k-1)*nx + i
                    ! ... Select seepage surface nodes
                    IF (ibc(ik,j)/10 == 3 .OR. ibc(ik,j)/10 == 5) THEN
                       IF(seeping(ik,j)) THEN
                          CALL gov_seep(gveqnm,gveqne,i,j,k)
                          ic = (k-1)*nx + i + (j-1)*nxz
                          ! ... Add in the capacitance terms. No source/sink term for seeping cell
                          gveqnm = gveqnm - (xm(ic)-xmoldt(ic))/delt
                          gveqne = gveqnm*h(ic)
                          ! ... Change sign so inflow is positive
                          gveqnm = -gveqnm
                          gveqne = -gveqne
                          aprint1(ic) = gveqnm
                          aprint2(ic) = gveqne
                          lprint1(ik,j) = 1
                          IF (ioptpr(2) == 0) THEN  ! ... 80 columns
                             WRITE(fubcf,9950) ltime, timedisply, i, j, k, &
                                  'Seepage', gveqnm, gveqne
                          ELSE                    ! ... 132 or 159 columns
                             WRITE(fubcf,9955) ltime, timedisply, i, j, k, &
                                  'Seepage', gveqnm, gveqne
                          END IF
                       END IF
                    END IF
                 END DO
              END DO
           END DO
           WRITE(fubcf2,2043)  &
                'Seepage Surface B.C. Flow Rates (at end of time step)',  &
                '(positive is into the region)'
2043       FORMAT(//tr40,a/tr50,a)
           WRITE(fubcf2,2042) 'Fluid Mass   (g/s)'
           CALL printar(2,aprint1,lprint1,fubcf2,1._kdp,24,000)
           WRITE(fubcf2,2042) 'Associated Advective Heat   (erg/s)'
           CALL printar(2,aprint2,lprint1,fubcf2,1._kdp,24,000)
        END IF
        IF(nhcond > 0) THEN        ! ... Heat conduction cells along base of region
           aprint1 = 0._kdp
           lprint1 = -1
           k = 1
           DO  j = 1,ny
              DO  i = 1,nx
                 ij = (j-1)*nx + i
                 ik = (k-1)*nx + i
                 ic = (j-1)*nxz + (k-1)*nx + i              
                 IF (MOD(ibc(ik,j),10) == 2) THEN
                    IF (irad) THEN
                       qcdtndum = cdtn(ij)*pi*(rsq(i+1)-rsq(i))
                    ELSE
                       qcdtndum = cdtn(ij)*dx(i)*dy(j)
                    END IF
                    aprint1(ic) = qcdtndum
                    lprint1(ik,j) = 1
                    IF (ioptpr(2) == 0) THEN
                       WRITE(fubcf,9050) ltime, timedisply, i, j, k, 'BsCnd',  &
                            qcdtndum, qcdtndum
                    ELSE
                       !                 132 or 159 columns
                       WRITE(fubcf,9055) ltime, timedisply, i, j, k, 'BslCond',  &
                            qcdtndum, qcdtndum
                    END IF
                 END IF
              END DO
           END DO
           WRITE(fubcf2,2043)  &
                'Heat Flux B.C. Flow Rates (at end of time step)',  &
                '(positive is into the region)'
           WRITE(fubcf2,2042) 'Heat   (erg/s)'
           CALL printar(2,aprint1,lprint1,fubcf2,1._kdp,24,000)
        END IF
        IF(nprecip > 0) THEN        ! ... Precipitation flux cells along top of region
           aprint1 = 0._kdp
           aprint2 = 0._kdp
           lprint1 = -1
           k = nz
           DO  j = 1,ny
              DO  i = 1,nx
                 ij = (j-1)*nx + i
                 ikp = (ktop(ij)-1)*nx + i
                 icp = (j-1)*nxz + (ktop(ij)-1)*nx + i
                 ! ... Select precipitation recharge cells or non-seeping recharge cells
                 IF (ibc(ikp,j) == 20 .OR. (ibc(ikp,j) == 50 .AND. .NOT.seeping(ikp,j))) THEN
                    IF (irad) THEN
                       qprecipdum = denflux(ij)*qprecip(ij)*pi*(rsq(i+1)-rsq(i))
                    ELSE
                       qprecipdum = denflux(ij)*qprecip(ij)*dx(i)*dy(j)
                    END IF
                    qhdum = ehflux(ij)*qprecipdum
                    aprint1(icp) = qprecipdum
                    aprint2(icp) = qhdum
                    lprint1(ikp,j) = 1
                    IF (ioptpr(2) == 0) THEN  ! ... 80 columns
                       WRITE(fubcf,9950) ltime, timedisply, i, j, k,  &
                            'Preci', qprecipdum, qhdum
                    ELSE                    ! ... 132 or 159 columns
                       WRITE(fubcf,9955) ltime, timedisply, i, j, k,  &
                            'Precip ', qprecipdum, qhdum
                    END IF
                 END IF
              END DO
           END DO
           WRITE(fubcf2,2043)  &
                'Precipitation Flux B.C. Flow Rates (at end of time step)',  &
                '(positive is into the region)'
           WRITE(fubcf2,2042) 'Fluid Mass   (g/s)'
           CALL printar(2,aprint1,lprint1,fubcf2,1._kdp,24,000)
           WRITE(fubcf2,2042) 'Associated Advective Heat   (erg/s)'
           CALL printar(2,aprint2,lprint1,fubcf2,1._kdp,24,000)
        END IF
     END IF
  END IF
  IF(ABS(print_dimno) > 0._kdp) THEN
     CALL print_control(print_dimno,time,ltime,tchg,timprdimno,lprint)
     IF(flag == 1 .OR. flag == 2) lprint = .FALSE.      ! ... no dimensionless nos. if error exit
     IF(lprint) THEN
        ! ... Compute and print Peclet and Nusselt numbers for all interior active cells
        CALL dimnum
        ! ... Set mask for Peclet to only print liquid water phase values
        lprint1 = -1
        DO  j = 1,ny
           DO  k = nz,1,-1
              DO  i = 1,nx
                 ik = (k-1)*nx + i
                 ic = (j-1)*nxz + (k-1)*nx + i
                 IF(ibc(ik,j) /= -1 .AND. ind(ic) /=2 .AND. ind(ic) /= 3 .AND. ind(ic) /= 5) THEN
                    lprint1(ik,j) = 1
                 END IF
              END DO
           END DO
        END DO
        IF(ltime == 0) THEN
           IF(.NOT.irad) THEN
              WRITE(fudimno,2005) rxlbl//'-Direction Node Coordinates  '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fudimno,unitfac(12),12,nx)
              WRITE(fudimno,2005) 'Y-Direction Node Coordinates  '//TRIM(unitlabel(13))
              CALL printar(1,yy,ibc,fudimno,unitfac(13),12,ny)
           ELSE IF(irad) THEN          ! ... Cylindrical aquifer geometry
              WRITE(fudimno,2005) rxlbl//'-Direction Node Coordinate '//TRIM(unitlabel(12))
              CALL printar(1,xx,ibc,fudimno,unitfac(12),12,nr)
           END IF
           WRITE(fudimno,2005) 'Z-Direction Node Coordinates   '//TRIM(unitlabel(14))
           CALL printar(1,zz,ibc,fudimno,unitfac(14),12,nz)
        END IF
        IF (iyr == 2) THEN
           WRITE(fudimno,9015) 'Time Step No. ',ltime,'; Simulation time: ',year(time),' (yr)'
        ELSE
           WRITE(fudimno,9015) 'Time Step No. ',ltime,'; Simulation time: ',time,' (s)'
        END IF
        IF (ioptprta(27) == 1) THEN
           WRITE(fudimno,9010) 'Cell Thermal Peclet Number','(-)'
           CALL printar(2,pe,lprint1,fudimno,1._kdp,24,000)
        END IF
        IF (ioptprta(28) == 1) THEN
           WRITE(fudimno,9010) 'Cell Nusselt Number','(-)'
           CALL printar(2,nu,ibc,fudimno,1._kdp,24,000)
        END IF
     END IF
  END IF
  ! ... Set the next time for printout
  prtimenxt = MIN(tchg, timprp, timprh, timprt, timprsat, timprd, timprv, timprpot,  &
       timprvel, timprbcf, timprsrc, timprprop, timprpor, timprperm, timprbal, timprdump,  &
       timprresid, timprdimno, timprpscal, timprpvec, timprtemp)
  DEALLOCATE(lprint1, hvrtavg, pvrtavg, swvrtavg, tcvrtavg, dumarry, aprint1, aprint2)
9005 FORMAT (1P, 8G10.3/(8G10.3))
9010 FORMAT (/tr10, '--- ', a, ' ---  ', a)
9015 FORMAT (/tr5,a,i5,a,1pg14.7,a)
9020 FORMAT (tr5,'Max magnitude residual ... ',1pg12.5,' at I,J,K ',3I4,  &
       '  N-R criterion <',1pg10.3)
9035 FORMAT (tr5,'Maximum value ... ',1pg11.4)
9050 FORMAT (i4, 1P, 1X, g8.2, i4, i3, i4, 1X, a5, 2(10X,g10.2))
9055 FORMAT (i5, 1P, 2X, g9.3, 3I4, 2X, a7, 2(14X,g10.3))
9950 FORMAT (i4, 1P, 1X, g8.2, i4, i3, i4, 1X, a5, g10.2, tr10, g10.2)
9955 FORMAT (i5, 1P, 2X, g9.3, 3I4, 2X, a7, tr2, g10.3, tr14, g10.3)
END SUBROUTINE pdata
