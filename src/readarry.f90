SUBROUTINE readarry(narray)
  ! ... Purpose:  To read the input arrays by free or specified format
  ! ...     and to calculate values for arrays with "CALCULATE" specified
  ! ... Note:  Pressure and Enthalpy arrays are calculated in
  ! ...            another subroutine called from GDATA (if necessary)
  USE machine_constants, ONLY: kdp
  USE f_units
  USE bc
  USE parameters          !!$*** for flint
  USE control
  USE mesh
!!$  USE parameters
  USE units, ONLY: unitfac
  USE variables
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: narray      ! number of arrays to be read
  !
  TYPE :: f_depth
     INTEGER :: zone_no, fn_no
     REAL(kind=kdp) :: par1, par2, par3, par4
  END TYPE f_depth
!!$  CHARACTER(LEN=76), EXTERNAL :: uppercase
  CHARACTER (LEN=3) :: utop
  CHARACTER (LEN=4) :: array, &
       zper='ZPER', xspa='XSPA', yspa='YSPA', zspa='ZSPA',  &
       xper='XPER', yper='YPER', ther='THER', cond='COND',  &
       enth='ENTH', temp='TEMP', pres='PRES', spec='SPEC',  &
       rden='RDEN', comp='COMP', precip='PREC', tfluxn='TFLU'
  CHARACTER(LEN=80) :: fmt
  CHARACTER(LEN=76) :: dumstring
  CHARACTER(LEN=80) :: duminput
  INTEGER :: a_err, da_err, i, ic, icalcopt, icup, icxf, icxl, ij, ik,  &
       ifunct, ioptperm, ioptporo, ioptpres, iopttemp, irx, iz, j, k, kdum, kk, kodd, lll, ncols
  INTEGER :: dep_sca_opt, intyp, kpar, nzppts, nztpts, rd_err, zone_no
  REAL(KIND=kdp) :: a0, a1, alf, depth, depth1, depth2, dum1, flxleft, flxright, gam,  &
       lcnv,  &
       perm1, perm2, permtemp, poro1, poro2, porotemp, presstemp, ptopalt,  &
       permbot, permtop, porobot, porotop, tempbot, temptop,  &
       tempbotalt, temptemp, temp1, temp2, temptopalt, tflxleft, tflxright, xkmax, xratio, z_datum,  &
       zele1, zele2
  REAL(KIND=kdp) :: interp1
  REAL(KIND=kdp), DIMENSION(6) :: zt(6), tvz(6), zp(6), pvz(6)
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: dumxy
  INTEGER :: alloc_stat
  TYPE (f_depth), DIMENSION(:), ALLOCATABLE :: zone_prop
  INTERFACE
     FUNCTION uppercase(string) RESULT(outstring)
       CHARACTER(LEN=*), INTENT(IN) :: string
       CHARACTER(LEN=LEN(string)) :: outstring
     END FUNCTION uppercase
     SUBROUTINE rdspace(dumarray, array, fmt, utop, nl, kodd, convfac)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: dumarray
       CHARACTER(LEN=4), INTENT(IN) :: array
       CHARACTER(LEN=80), INTENT(IN) :: fmt
       CHARACTER(LEN=3), INTENT(IN) :: utop
       INTEGER, INTENT(IN) :: nl
       INTEGER, INTENT(OUT) :: kodd
       REAL(KIND=kdp), INTENT(IN) :: convfac
     END SUBROUTINE rdspace
     SUBROUTINE rd(dumarray, fmt, utop, kodparm, kodd, dumarray2)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: dumarray
       CHARACTER(LEN=80), INTENT(IN) :: fmt
       CHARACTER(LEN=3), INTENT(IN) :: utop
       INTEGER, INTENT(IN) :: kodparm
       INTEGER, INTENT(OUT) :: kodd
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT), OPTIONAL :: dumarray2
     END SUBROUTINE rd
     SUBROUTINE presboil(p,h,tc,ii,jj)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: p, h, tc
       INTEGER, INTENT(IN) :: ii, jj
     END SUBROUTINE presboil
  END INTERFACE
  !     ------------------------------------------------------------------
  ! ...   Cdtn - Conductive Heat Flux along the base of region
  ! ...   Dx -  cell sizes in X direction
  ! ...   Dy -  cell sizes in Y direction
  ! ...   Dz -  cell sizes in Z direction
  ! ...   H -  enthalpy
  ! ...   kodparm - parameter index
  ! ...    1 -  Porosity
  ! ...    2 -  X - Permeability
  ! ...    3 -  Y - Permeability
  ! ...    4 -  Z - Permeability
  ! ...    5 -  Thermal Conductivity
  ! ...    6 -  Specific Heat of Rock
  ! ...    7 -  Rock Density
  ! ...    8 -  Rock compressibility
  ! ...    9 -  Conductive Heat Flux along the Base of region
  ! ...    10 -  Enthalpy or Temperature
  ! ...    11 -  Pressure
  ! ...    12 -  X - spacing
  ! ...    13 -  Y - spacing
  ! ...    14 -  Z - spacing
  ! ...    18 -  Precipitation Flux along the top of region
  ! ...    19 -  Associated temperature for precipitation flux along the top of region
  ! ...   Kod(kodparm) =1 - array input as variable value array
  ! ...             =2 - array input as constant value array
  ! ...             =3 - array input as values constant for each row
  ! ...             =4 - array input as values constant for each column
  ! ...             =5 - array input as rocktype values
  ! ...             6 - array input as rocktype with function of temperature
  ! ...             7 - array input as rocktype with function of time
  ! ...             8 - array input as rocktype with function of depth and temperature
  ! ...         =50 -  code for writing Explorer plotfiles
  ! ...         =60 -  code for hydrostatic Pressure gradient
  ! ...         =70 -  code indicating P & H of all nodes are on the sat'd
  ! ...                  water or steam curve (hydrostatic P gradient for boiling
  ! ...                  point density)
  ! ...         =75 -  code indicating P & H of some columns are on the
  ! ...                  sat'd water or steam curve (hydrostatic P gradient for
  ! ...                  boiling point density)
  ! ...         =77 -  code indicating P & H of a node is on the sat'd
  ! ...                 water or steam curve (Note: no hydrostatic P gradient
  ! ...                 for boiling is implied)
  ! ...         =90 -  code indicating Y and/or Z permeabilities
  ! ...                 are identical to X permeabilities
  ! ...   Kodt  =0 if enthalpy data are input for all nodes
  ! ...          4 if enthalpy data are input for only a subset of nodes
  ! ...          1 if temperature data are input for all nodes
  ! ...          2 if temperature data are input for only a subset of nodes
  ! ...          3 if temperature data are input with boiling point options
  ! ...   P -  pressure
  ! ...   Phi - porosity
  ! ...   Ptop - Pressure at top of domain -used for auto-P opt otherwise =-99
  ! ...   Xk -  X permeability
  ! ...   Xkc - Thermal conductivity
  ! ...   Yk -  Y permeability
  ! ...   Ykfactor - value read in to compute Y permeabilites from X perm
  ! ...   Zk -  Z permeability
  ! ...   Zkfactor - value read in to compute Z permeabilites from X perm
  !
  ! ... Local Variables:
  ! ...   alf - variable used in calculating radial coordinates
  ! ...   flxleft - heat conduction calculation variable
  ! ...   flxright - heat conduction calculation variable
  ! ...   fmt - format
  ! ...   ioptpres - pressure calculation variable
  ! ...   iopttemp - temperature calculation variable
  ! ...   tempbot - temperature calculation variable
  ! ...   temptemp - temperature calculation variable
  ! ...   temptop - temperature calculation variable
  ! ...   xratio - heat conduction calculation variable
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.3 $//$Date: 2007/10/16 21:17:22 $'
  !     ------------------------------------------------------------------
  ALLOCATE(dumxy(nxy), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: readarray'
    ierr(180) = .TRUE.
    RETURN
  END IF
  loop310:  DO  kk=1,narray
     ! ...   read in character string containing array name and strip leading blanks
     READ(fuins,'(a)') duminput
     duminput = ADJUSTL(duminput)
     ! ... split into name of array variable and units string
     array = duminput(1:4)
     dumstring = duminput(5:80)
     WRITE(fuclog,'(10X,2A)') 'Reading Input Array - ', array
     WRITE(fupdef,'(10X,2A)') 'Reading Input Array - ', array
     ! ... determine units (SI, cgs, etc.) for this array variable
     array = uppercase(array)
     CALL getunits(array, dumstring, kodd)
     ! ... read in character string containing format and strip leading blanks
     READ(fuins,'(a)') fmt
     fmt = ADJUSTL(fmt)
     duminput = uppercase(fmt(1:10))
     ! ...  specified format may be input or one of:
     ! ...         'FREE' = free format
     ! ...         'CONS' = constant value
     ! ...         'ROWC' = variables are constant for each row
     ! ...         'COLC' = variables are constant for each column
     ! ...         'NODE' = variables are input node by node
     ! ...         'CALC' = calculate values
     ! ...         'ROCK' = variables are input by rocktype
     ! ...         'RTFDEPTH' = calculate values as f(depth) for variables input by rocktype
     ! ...         'RTFTEMP' = calculate values as f(Temperature) for variables input by rocktype
     ! ...         'RTFTIME' = calculate values as f(time) for variables input by rocktype
     !
     ! ...    if input version is VER2 or VER3 and
     ! ...      if fmt = 'FREE' or 'ROWC'
     ! ...      then 'TOP' or 'BOTTOM' required for
     ! ...      the following arrays: ZSPA,PHI,XK,YK,ZK,XKC,PHFWT,DF,BETA,H,P
     utop = '   '
     IF ((informt >= 2) .AND. (duminput(1:4) == 'FREE' .OR. duminput(1:4) == 'ROWC') .AND.  &
          (array == zspa .OR. array == 'PORO' .OR. array == xper .OR.   &
          array == yper .OR. array == zper .OR. array == ther .OR.   &
          array == spec .OR. array == rden .OR. array == comp .OR.   &
          array == enth .OR. array == temp .OR. array == pres)) THEN
        ! ...       read in character string containing 'TOP' or 'BOT' and
        ! ...         strip leading blanks
        READ(fuins,'(a)') dumstring
        dumstring = ADJUSTL(dumstring)
        utop = uppercase(dumstring(1:3))
        IF (utop /= 'TOP' .AND. utop /= 'BOT') THEN
           WRITE (fustdout,9050) array, fmt(1:7), array, fmt(1:7)
           WRITE (fuclog,9050) array, fmt(1:7), array, fmt(1:7)
9050       FORMAT (/, '**** INPUT ERROR: Reading array ', a, ' *****', /, 5X,  &
                a, ' format for ', a, ' requires keyword TOP or BOTTOM on line following ', a)
           ierr(34) = .TRUE.
        END IF
     ELSE
        utop = '   '
     END IF
     ! ... Read keyword-data block if the format specifier is neither CALC nor RTFDEPTH
     IF (duminput(1:4) /= 'CALC' .AND. duminput(1:8) /= 'RTFDEPTH') THEN
        IF (array == xspa) THEN        ! ...                 ----- X Spacing -----
           CALL rdspace(dx, array, fmt, utop, nx, kod(12), unitfac(12))
        ELSE IF (array == yspa) THEN   ! ...                 ----- Y Spacing -----
           IF(.NOT.irad .OR. informt <= 2) THEN   !VER2 requires YPERM for cylindrical coord
              CALL rdspace(dy, array, fmt, utop, ny, kod(13), unitfac(13))
           ELSE
              WRITE (fustdout,9006) array
              WRITE (fuclog,9006) array
9006          FORMAT (10X, '**** INPUT ERROR: ', a4,  &
                   ' is not accepted array name with cylindrical coordinates***')
              ierr(46) = .TRUE.
           END IF
        ELSE IF (array == zspa) THEN   ! ...                 ----- Z Spacing -----
           CALL rdspace(dz, array, fmt, utop, nz, kod(14), unitfac(14))
        ELSE IF (array == 'PORO') THEN    ! ...       ----- Porosity -----
           CALL rd(phi, fmt, utop, 1, kod(1))
        ELSE IF (array == xper) THEN      ! ...       ----- X Permeability -----
           CALL rd(xk, fmt, utop, 2, kod(2))
        ELSE IF (array == yper) THEN      ! ...       ----- Y Permeability -----
            IF(.NOT.irad .OR. informt <= 2) THEN   !VER2 requires YPERM for cylindrical coord
               CALL rd(yk, fmt, utop, 3, kod(3))
            ELSE
               WRITE (fustdout,9106) array
               WRITE (fuclog,9106) array
 9106          FORMAT (10X, '**** INPUT ERROR: ', a4,  &
                    ' is not accepted array name with cylindrical coordinates***')
               ierr(46) = .TRUE.
            END IF
        ELSE IF (array == zper) THEN      ! ...       ----- Z Permeability -----
           CALL rd(zk, fmt, utop, 4, kod(4))
        ELSE IF (array == ther) THEN      ! ...       ----- Thermal Conductivity -----
           CALL rd(xkc, fmt, utop, 5, kod(5))
        ELSE IF (array == spec) THEN      ! ...       ----- Specific Heat (Heat Capacity) -----
           CALL rd(phfwt, fmt, utop, 6, kod(6))
        ELSE IF (array == rden) THEN      ! ...       ----- Rock Density -----
           CALL rd(df, fmt, utop, 7, kod(7))
        ELSE IF (array == comp) THEN      ! ...       ----- Rock Compressibility -----
           CALL rd(beta, fmt, utop, 8, kod(8))
        ELSE IF (array == cond) THEN      ! ...       ----- Conductive Heat Flux along Base -----
           IF (duminput(1:4) == 'FREE') THEN     ! ... free format
              READ(fuins,*) (cdtn(ij), ij=1, nxy)
              ! ...            convert from user specified units to internal cgs units
              DO  ij = 1,nxy
                 cdtn(ij) = cdtn(ij)*unitfac(9)
              END DO
              kod(9) = 1
              CYCLE loop310
           ELSE IF (duminput(1:4) == 'CONS') THEN     ! ... constant value
              READ(fuins,*) dum1
              ! ...            convert from user specified units to base units (cgs)
              ! ...              and assign value to array
              dum1 = dum1*unitfac(9)
              DO  ij = 1, nxy
                 cdtn(ij) = dum1
              END DO
              kod(9) = 2
              CYCLE loop310
           ELSE IF (duminput(1:4) == 'ROWC' .OR. duminput(1:4) == 'COLC' .OR.  &
                duminput(1:4) == 'ROCK') THEN
              ! ... constant value for each row or column or by rock type
              GO TO 990     ! ... not allowed
           ELSE IF (duminput(1:4) == 'NODE') THEN     ! ... individual node input
              kod(9) = 1
              DO  kdum = 1, nxyz + 1
                 READ(fuins,*) i, j, k, dum1
                 ! ...                     stop input if I,J,K=0
                 IF (i == 0 .AND. j == 0 .AND. k == 0) CYCLE loop310
                 IF (k /= 1) THEN
                    WRITE (fustdout,'(a)') '**** INPUT ERROR: K must be 1 for CONDUCTION ***'
                    WRITE (fuclog,'(/a)') '**** INPUT ERROR: K must be 1 for CONDUCTION ***'
                    ierr(27) = .TRUE.
                 END IF
                 ij = (j-1)*nx + i
                 ! ...           convert from user specified units to base units (cgs)
                 ! ...             and assign value to node
                 cdtn(ij) = dum1*unitfac(9)
              END DO
           ELSE
              ! ... specified format
              READ(fuins, fmt) (cdtn(ij), ij=1, nxy)
              ! ...            convert from user specified units to base units (cgs)
              DO  ij = 1, nxy
                 cdtn(ij) = cdtn(ij)*unitfac(9)
              END DO
              kod(9) = 1
           END IF
        ELSE IF (array == enth .OR. array == temp) THEN    ! ... ----- Enthalpy or Temperature ----
           IF (array == temp) THEN
              kodt = 1           ! ... temperature data are input for all nodes
              IF (duminput(1:4) == 'NODE') kodt = 2     ! ... temperature data are input for a subset of nodes
              CALL rd(tcassoc, fmt, utop, 15, kod(10), ehassoc)
           ELSE
              kodt = 0           ! ... enthalpy data are input for all nodes
              IF (duminput(1:4) == 'NODE') kodt = 4     ! ... enthalpy data are input for a subset of nodes
              CALL rd(ehassoc, fmt, utop, 10, kod(10))
           END IF
           IF (array == temp .AND. kodt >= 1) THEN           ! ... check that temperature is 
              !                                                    within limits of table
              DO  lll = 1, npicmax
                 ic = npic(lll)
                 IF (tcassoc(ic) > 1400.0_kdp .OR. tcassoc(ic) <= 0.0_kdp) THEN
                    CALL ictoijk(ic,i,j,k,nx,nz)
                    WRITE (fustdout,9030) 'T', ic, i, j, k, 'T', tcassoc(ic)
                    WRITE (fuclog,9030) 'T', ic, i, j, k, 'T', tcassoc(ic)
9030                FORMAT (' **** INPUT ERROR: ', a, ' value outside table limits at IC =', i4,  &
                         ';   I,J,K = ',3I4, ' *****'/tr11,a,' = ',1pg12.5)
                    ierr(30) = .TRUE.
                 END IF
              END DO
           ELSE         ! ... check that enthalpy is within limits of table
              DO  lll = 1, npicmax
                 ic = npic(lll)
                 IF (ehassoc(ic) > 52.0e9_kdp .OR. ehassoc(ic) < 0.01e9_kdp) THEN
                    CALL ictoijk(ic,i,j,k,nx,nz)
                    WRITE (fustdout,9030) 'H', ic, i, j, k, 'H', ehassoc(ic)
                    WRITE (fuclog,9030) 'H', ic, i, j, k, 'H', ehassoc(ic)
                    ierr(31) = .TRUE.
                 END IF
              END DO
           END IF
        ELSE IF (array == pres) THEN           ! ...    ----- Pressure -----
           ptop = -99._kdp
           CALL rd(p, fmt, utop, 11, kod(11))
           ! ...     check that pressure is within limits of table
           DO  lll = 1, npicmax
              ic = npic(lll)
              IF(unconfined) THEN
                 IF (p(ic) > 1.0e10_kdp) THEN
                    CALL ictoijk(ic,i,j,k,nx,nz)
                    WRITE (fustdout,9030) 'P', ic, i, j, k, 'P', p(ic)
                    WRITE (fuclog,9030) 'P', ic, i, j, k, 'P', p(ic)
                    ierr(47) = .TRUE.
                 END IF
              ELSE
                 IF (p(ic) > 1.0e10_kdp .OR. p(ic) < 0.5e6_kdp) THEN
                    CALL ictoijk(ic,i,j,k,nx,nz)
                    WRITE (fustdout,9030) 'P', ic, i, j, k, 'P', p(ic)
                    WRITE (fuclog,9030) 'P', ic, i, j, k, 'P', p(ic)
                    ierr(47) = .TRUE.
                 END IF
              END IF
           END DO
        ELSE IF (array == precip) THEN      ! ... ----- Precipitation Flux at top of region -----
           IF (duminput(1:4) == 'FREE') THEN     ! ... spatially variable distribution
              READ(fuins,*) (qprecip(ij), ij=1, nxy)
              ! ... convert from user specified units to internal units (cgs)
              DO  ij = 1,nxy
                 qprecip(ij) = qprecip(ij)*unitfac(18)
              END DO
              kod(18) = 1
           ELSE IF (duminput(1:4) == 'CONS') THEN     ! ... constant value
              READ(fuins,*) dum1
              ! ... convert from user specified units to internal units (cgs)
              ! ...       and assign value to array
              dum1 = dum1*unitfac(18)
              DO  ij = 1,nxy
                 qprecip(ij) = dum1
              END DO
              kod(18) = 2
           ELSE IF (duminput(1:4) == 'ROWC' .OR. duminput(1:4) == 'COLC' .OR.  &
                duminput(1:4) == 'ROCK') THEN
              ! ... constant value for each row or column or by rock type
              GO TO 990     ! ... This format not available
           ELSE IF (duminput(1:4) == 'NODE') THEN    ! ... individual node input
              kod(18) = 1
              DO  kdum = 1, nxyz + 1
                 READ(fuins,*) i, j, k, dum1
                 IF (i == 0 .AND. j == 0 .AND. k == 0) EXIT       ! ... stop input if I,J,K=0
                 IF (k /= nz) THEN
                    WRITE (fustdout,'(a)') '**** INPUT ERROR: K must be nz for precipitation ***'
                    WRITE (fuclog,'(/a)') '**** INPUT ERROR: K must be nz for precipitation ***'
                    ierr(48) = .TRUE.
                 END IF
                 ij = (j-1)*nx + i
                 ! ... convert from user specified units to internal units (cgs)
                 ! ...       and assign value to node
                 qprecip(ij) = dum1*unitfac(18)
              END DO
           ELSE              ! ... specified format
              READ(fuins,fmt) (qprecip(ij), ij=1, nxy)
              ! ... convert from user specified units to internal units (cgs)
              DO  ij = 1, nxy
                 qprecip(ij) = qprecip(ij)*unitfac(18)
              END DO
              kod(18) = 1
           END IF
        ELSE IF (array == tfluxn) THEN      ! ... Associated Temperature for Precipitation Flux
           IF (duminput(1:4) == 'FREE') THEN     ! ... spatially variable distribution
              READ(fuins,*) (tflux(ij), ij=1,nxy)
              ! ... Already in base units (Deg.C)
              kod(19) = 1
              GO TO 820
           ELSE IF (duminput(1:4) == 'CONS') THEN         ! ... constant value
              READ(fuins,*) dum1
              ! ... assign value to array
              DO  ij = 1, nxy
                 tflux(ij) = dum1
              END DO
              kod(19) = 2
              GO TO 820
           ELSE IF (duminput(1:4) == 'ROWC' .OR. duminput(1:4) == 'COLC' .OR.  &
                duminput(1:4) == 'ROCK') THEN
              ! ... constant value for each row or column
              GO TO 990     ! ... not allowed
           ELSE IF (duminput(1:4) == 'NODE') THEN         ! ... individual node input
              kod(19) = 1
33            READ(fuins,*) i, j, k, dum1
              IF (i == 0 .AND. j == 0 .AND. k == 0) GO TO 820      ! ... terminate input
              IF (k /= nz) THEN
                 WRITE (fustdout,'(a)') '**** INPUT ERROR: K must be nz for TFLUX ***'
                 WRITE (fuclog,'(/a)') '**** INPUT ERROR: K must be nz for TFLUX ***'
                 ierr(28) = .TRUE.
              END IF
              ij = (j-1)*nx + i
              ! ... already in internal units (cgs), assign value to node
              tflux(ij) = dum1
              GO TO 33
           ELSE              ! ... specified format
              READ(fuins,fmt) (tflux(ij), ij=1,nxy)
              ! ... already in internal units (cgs)
              kod(19) = 1
           END IF
820        CONTINUE
           ! ...     check that Temperature is within limits of table
           DO  ij=1,nxy
              IF (tflux(ij) > 1400._kdp .OR. tflux(ij) <= 0._kdp) THEN
                 j = ij/nx + 1
                 i = ij - (j-1)*nx
                 WRITE (fustdout,9031) 'T', ij, i, j, 'Tflux', tflux(ij)
                 WRITE (fuclog,9031) 'T', ij, i, j, 'Tflux', tflux(ij)
9031             FORMAT ('**** INPUT ERROR: ', a,' outside table limits at IJ =', i4,  &
                      ';   I,J = ', 2I4, ' *****' /tr11,a,' = ',1pg12.5)
                 ierr(32) = .TRUE.
              END IF
           END DO
        ELSE
           ! ...   #######  array name is not valid - write error message #######
           WRITE (fustdout,9005) array
           WRITE (fuclog,9005) array
9005       FORMAT (tr10, '**** INPUT ERROR: ',a4,' is invalid array name ***')
           ierr(49) = .TRUE.
        END IF
     ELSE IF (duminput(1:4) == 'CALC') THEN        ! ...   Calculate Values for Arrays
        IF (array == xspa .AND. irad) THEN
           READ(fuins,*) xx(1)
           ! ...          XX(1) - radius of first node from origin
           ! ...          RSQ(1) - well radius
           ! ...          RSQ(NXX) - max reservoir radius
           ! ...      convert user specified units to internal base units (cm)
           xx(1) = xx(1)*unitfac(12)
           alf = (rsq(nxx)/xx(1))**(1.0_kdp/(DBLE(nx)-0.5_kdp))
           DO  i = 2, nx
              xx(i) = xx(i-1)*alf
           END DO
           ! ...      assign dummy value to node spacing - these are calculated in gdata
           dx(1) = -1.0_kdp
           kod(12) = 1
! ***  deactivate the permeability, superceded by rock type function of depth ***
!!$        ELSE IF (array == xper) THEN           ! ... ----- X Permeability -----
!!$           READ(fuins,*) ioptperm
!!$           ! ... Permeability option 1: linear gradient from top to bottom
!!$           ! ...          for permeability option 1 set KOD(2)=3 for constant row values
!!$           ! ...             except where NX=1  set KOD(2)=1
!!$           IF (ioptperm == 1) THEN
!!$              READ(fuins,*) permtop, permbot
!!$              permtop = permtop*unitfac(2)
!!$              permbot = permbot*unitfac(2)
!!$              kod(2) = 3
!!$              IF (nx == 1) kod(2) = 1
!!$              DO  k = 1, nz
!!$                 permtemp = permtop + (zz(nz)-zz(k))*(permbot-permtop)/(zz(nz)-zz(1))
!!$                 DO  j = 1, ny
!!$                    icxf = (k-1)*nx + 1 + (j-1)*nxz
!!$                    icxl = icxf + nx - 1
!!$                    DO  ic = icxf, icxl
!!$                       xk(ic) = permtemp
!!$                    END DO
!!$                 END DO
!!$              END DO
!!$           ELSE IF (ioptperm == 2) THEN
!!$              READ(fuins,*) permtop, permbot
!!$              ! ... Permeability option 2: exponential decrease of permeability with depth
!!$              ! ...      based on power of 10
!!$              ! ...      set KOD(2)=3 for constant row values
!!$              ! ...           except if NX=1, set KOD(2)=1
!!$              permtop = permtop*unitfac(2)
!!$              permbot = permbot*unitfac(2)
!!$              gam = LOG10(permtop/permbot)/(zz(nz)-zz(1))
!!$              kod(2) = 3
!!$              IF (nx == 1) kod(2) = 1
!!$              DO  k = 1,nz
!!$                 permtemp = permtop*10**(gam*zz(k))
!!$                 DO  j = 1,ny
!!$                    icxf = (k-1)*nx + 1 + (j-1)*nx*nz
!!$                    icxl = icxf + nx - 1
!!$                    DO  ic =icxf,icxl
!!$                       xk(ic) = permtemp
!!$                    END DO
!!$                 END DO
!!$              END DO
!!$           ELSE IF (ioptperm == 3) THEN
!!$              READ(fuins,*) permtop, permbot
!!$              ! ... Permeability option 3: logarithmic variation of permeability with log of depth
!!$              ! ...     set KOD(2)=1 for variable values
!!$              ! ... Depth is from land surface 
!!$              depth1 = zls(1)-zz(nz)
!!$              depth2 = zls(1)-zz(1)
!!$              ! ... Convert from user specified units to base units (cgs)
!!$              permtop = permtop*unitfac(2)
!!$              permbot = permbot*unitfac(2)
!!$              a0 = (LOG10(depth1)*LOG10(permbot)-LOG10(depth2)*LOG10(permtop))/  &
!!$                                  LOG10(depth1/depth2)
!!$              a1 = LOG10(permtop/permbot)/LOG10(depth1/depth2)
!!$              kod(2) = 1
!!$              DO  j=1,ny
!!$                 DO i=1,nx
!!$                    ij = (j-1)*nx + i
!!$                    DO  k=1,nz
!!$                       ic = (k-1)*nx + i + (j-1)*nxz
!!$                       depth = (zls(ij)-zz(k))         ! ... (cm)
!!$                       IF(depth >= 0._kdp) THEN
!!$                          xk(ic) = 10**(a0+a1*LOG10(depth))
!!$                       ELSE
!!$                          xk(ic) = 0._kdp
!!$                       END IF
!!$                    END DO
!!$                 END DO
!!$              END DO
!!$           ELSE IF (ioptperm == 4) THEN
!!$              READ(fuins,*) permtop, permbot
!!$              ! ... Permeability option 4: logarithmic variation of permeability with depth
!!$              ! ...    from Williams and Narasimhan (1989)
!!$              ! ...          set KOD(2)=1 for variable values
!!$              ! ... Depth is from land surface 
!!$              depth1 = zls(1)-zz(nz)
!!$              depth2 = zls(1)-zz(1)
!!$              ! ... Convert from user specified units to base units (cgs)
!!$              permtop = permtop*unitfac(2)
!!$              permbot = permbot*unitfac(2)
!!$              a0 = (depth1*LOG10(permbot)-depth2*LOG10(permtop))/(depth1-depth2)
!!$              a1 = LOG10(perm1/perm2)/(depth1-depth2)
!!$              kod(2) = 1
!!$              DO  j=1,ny
!!$                 DO i=1,nx
!!$                    ij = (j-1)*nx + i
!!$                    DO  k=1,nz
!!$                       ic = (k-1)*nx + i + (j-1)*nxz
!!$                       depth = (zls(ij)-zz(k))          ! ... (cm)
!!$                       xk(ic) = 10**(a0+a1*depth)
!!$                    END DO
!!$                 END DO
!!$              END DO
!!$           ELSE
!!$              GO TO 990
!!$           END IF
           ! ...                 ----- Y Permeability -----
           ! ...             compute using YKFACTOR*X-permeability
        ELSE IF (array == yper) THEN
           READ(fuins,*) ykfactor
           kod(3) = 90
           DO  ic = 1, nxyz
              yk(ic) = ykfactor*xk(ic)
           END DO
           ! ...              set the Y-Perm values for Rock types
           DO  lll = 1, nrxtype
              irx = irxused(lll)
              rxparm(3,irx) = ykfactor*rxparm(2, irx)
              ykrxkxfac(irx) = ykfactor
           END DO
           ! ...                 ----- Z Permeability -----
           ! ...             compute using ZKFACTOR*X-permeability
        ELSE IF (array == zper) THEN
           READ(fuins,*) zkfactor
           kod(4) = 90
           DO  ic = 1, nxyz
              zk(ic) = zkfactor*xk(ic)
           END DO
           ! ...              set the Z-Perm values for Rock types
           DO  lll = 1, nrxtype
              irx = irxused(lll)
              rxparm(4,irx) = zkfactor*rxparm(2,irx)
              zkrxkxfac(irx) = zkfactor
           END DO
! ***  deactivate the permeability, superceded by rock type function of depth ***
!!$        ELSE IF (array == 'PORO') THEN          ! ... ----- Porosity -----
!!$           READ(fuins,*) ioptporo, porotop, porobot
!!$           ! ... Porosity option 1: linear gradient from top to bottom
!!$           ! ...          for porosity option 1 set KOD(1)=3 for constant row values
!!$           ! ...             except where NX=1  set KOD(1)=1
!!$           IF (ioptporo == 1) THEN
!!$              kod(1) = 3
!!$              IF (nx == 1) kod(1) = 1
!!$              DO  k = 1, nz
!!$                 porotemp = porotop + (zz(nz)-zz(k))*(porobot-porotop)  &
!!$                      /(zz(nz)-zz(1))
!!$                 DO  j = 1, ny
!!$                    icxf = (k-1)*nx + 1 + (j-1)*nx*nz
!!$                    icxl = icxf + nx - 1
!!$                    DO  ic = icxf, icxl
!!$                       phi(ic) = porotemp
!!$                    END DO
!!$                 END DO
!!$              END DO
!!$              ! ... Porosity option 2: exponential decrease of porosity with depth
!!$              ! ...          for porosity option 2 set KOD(10)=3 for constant row values
!!$              ! ...                 ----- Porosity -----
!!$              ! ...             except where NX=1  set KOD(10)=1
!!$           ELSE IF (ioptporo == 2) THEN
!!$              gam = -LOG(porobot/porotop)
!!$              kod(1) = 3
!!$              IF (nx == 1) kod(1) = 1
!!$              DO  k = 1, nz
!!$                 porotemp = porotop* EXP(-gam*(zz(nz)-zz(k))/(zz(nz)-zz(1)))
!!$                 DO  j = 1, ny
!!$                    icxf = (k-1)*nx + 1 + (j-1)*nx*nz
!!$                    icxl = icxf + nx - 1
!!$                    DO  ic = icxf, icxl
!!$                       phi(ic) = porotemp
!!$                    END DO
!!$                 END DO
!!$              END DO
!!$           ELSE
!!$              GO TO 990
!!$           END IF
        ELSE IF (array == cond) THEN
           ! ...                 ----- Conductive Heat Flux along the Base -----
           ! ...     compute Conductive Heat input using a linear interpolation in the
           ! ...       X-direction.  In WELLALLO the flux (mW/m^2)=(erg/s-cm^2) is
           ! ...       multiplied by the
           ! ...       area (cm^2) of each cell to give a rate (erg/s)
           READ(fuins,*) icalcopt, flxleft, flxright
           ! ...         option 1 (presently the only option)
           IF (icalcopt == 1) THEN
              ! ...                if two values are equal
              IF (flxleft == flxright) THEN
                 kod(9) = 2
                 DO  ij = 1, nxy
                    cdtn(ij) = flxleft
                 END DO
              ELSE
                 ! ...          check for input error using variable heat flux with NX=1
                 ! ...             this will cause division by 0
                 IF (nx == 1) THEN
                    WRITE (fustdout, '(/,2A)')  &
                         '**** INPUT ERROR: Spatially variable conductive heat flux is not ',  &
                         'appropriate for NX=1 ***'
                    WRITE (fuclog, '(/,2A)')  &
                         '**** INPUT ERROR: Spatially variable conductive heat flux is not ',  &
                         'appropriate for NX=1 ***'
                    ierr(50) = .TRUE.
                    !..                  return
                 END IF
                 kod(9) = 1
                 DO  i = 1, nx
                    xratio = (xx(i)-xx(1))/(xx(nx)-xx(1))
                    DO  j = 1, ny
                       ij = (j-1)*nx + i
                       cdtn(ij) = ((1.0_KDP-xratio)*flxleft+xratio*flxright)
                    END DO
                 END DO
              END IF
              ! ...            convert from user specified units to base units (cgs)
              DO  ij = 1, nxy
                 cdtn(ij) = cdtn(ij)*unitfac(9)
              END DO
              CYCLE loop310
           END IF
           GO TO 990
        ELSE IF (array == temp) THEN        ! ...                 ----- Temperature -----
           ! ...      NOTE:  no unit conversion is necessary for Temperature since
           ! ...             deg C is the only choice at this time
           READ(fuins,*) iopttemp, temptop, tempbot
           IF (iopttemp == 1) THEN
              ! ... Temperature option 1: linear gradient from top temperature TEMPTOP to 
              ! ...     bottom temperature TEMPBOT
              kodt = 1
              kod(10) = 3           ! ... constant row values 
              IF (nx == 1) kod(10) = 1
              DO  k = 1,nz
                 temptemp = temptop + (zz(nz)-zz(k))*(tempbot-temptop)/(zz(nz)-zz(1))
                 DO  j = 1, ny
                    icxf = (k-1)*nx + 1 + (j-1)*nx*nz
                    icxl = icxf + nx - 1
                    DO  ic = icxf, icxl
                       tcassoc(ic) = temptemp
                    END DO
                 END DO
              END DO
           ELSE IF (iopttemp == 11) THEN
              ! ... Temperature option 11: piecewise linear temperature profile
              ! ...      from top temperature TEMPTOP to bottom temperature TEMPBOT 
              ! ...      with up to 5 linear segments
              kodt = 1
              kod(10) = 3              ! ... constant row values
              IF (nx == 1) kod(10) = 1
              ! ... z values are input from top to bottom as depth values and are 
              ! ...      indexed to relate to k increasing from
              ! ...           bottom to top of region
              READ(fuins,*) nztpts,(zt(i),tvz(i),i=nztpts,1,-1)
              IF(nztpts > 6) THEN
                 ierr(67) = .TRUE.
              END IF
              ! ...      convert from depth to z-coordinate and into cm
              DO  i=1,nztpts
                 zt(i) = zz(nz)+0.5_kdp*dz(nz) - zt(i)*unitfac(14)
              END DO
              DO  k=1,nz
                 temptemp = interp1(zz(k),nztpts,zt,tvz)
                 DO  j=1,ny
                    DO i=1,nx
                       ik = (k-1)*nx + i
                       ic = (k-1)*nx + i + (j-1)*nxz
                       ! ... Only load the active cells
                       IF (np(ik,j) == 0) THEN
                          tcassoc(ic) = 0._kdp
                       ELSE
                          tcassoc(ic) = temptemp
                       END IF
                    END DO
                 END DO
              END DO
           ELSE IF (iopttemp == 2) THEN
              ! ... Temperature option 2: linear gradient from top temperature,
              ! ...     TEMPTOP, calculated using input gradient value in TEMPBOT
              kodt = 1
              kod(10) = 1
              DO  j = 1, ny
                 DO  i = 1, nx
                    icup = 0
                    DO  k = nz, 1, -1
                       ic = (k-1)*nx + i + (j-1)*nx*nz
                       ik = (k-1)*nx + i
                       ! ... skip past nodes not in the domain
                       IF (np(ik,j) == 0) THEN
                          tcassoc(ic) = 0._kdp
                          CYCLE
                       END IF
                       ! ... assign input temperature to top cell in domain in each column
                       IF (icup == 0) THEN
                          tcassoc(ic) = temptop
                       ELSE
                          ! ... calc temperature at lower nodes using input thermal gradient
                          tcassoc(ic) = tcassoc(icup) + tempbot*0.5_kdp*(dz(k+1)+dz(k))/1.e5_kdp
                       END IF
                       icup = ic
                    END DO
                 END DO
              END DO
           ELSE IF (iopttemp == 21) THEN
              ! ... Temperature option 21: Similar to option 2 except that the
              ! ...          upper temperature values are not uniform, although the
              ! ...          linear gradient is uniform.  Temperature values must be
              ! ...          input for all columns in mesh.  This value is assigned to the
              ! ...          uppermost active node
              kodt = 1
              kod(10) = 1
              READ(fuins,*) (dumxy(ij), ij=1, nxy)
              DO  j = 1, ny
                 DO  i = 1, nx
                    ij = (j-1)*nx + i
                    icup = 0
                    DO  k = nz, 1, -1
                       ic = (k-1)*nx + i + (j-1)*nx*nz
                       ik = (k-1)*nx + i
                       ! ...             skip past nodes not in the domain
                       IF (np(ik,j) == 0) THEN
                          tcassoc(ic) = 0._kdp
                          CYCLE
                       END IF
                       ! ... assign input temperature to top cell in domain in each column
                       IF (icup == 0) THEN
                          tcassoc(ic) = dumxy(ij)
                       ELSE
                          ! ... calc temperature at lower nodes using input thermal gradient
                          tcassoc(ic) = tcassoc(icup) + tempbot*0.5_kdp*(dz(k+1)+dz(k))/1.e5_kdp
                       END IF
                       icup = ic
                    END DO
                 END DO
              END DO
           ELSE IF (iopttemp == 22) THEN
              ! ... Temperature option 22: Similar to option 2 except that the
              ! ...          linear gradient is only for specified columns (I,J)
              ! ...          the gradient is from the top temperature,
              ! ...          TEMPTOP, calculated with input gradient value, TEMPBOT
! *** This option is currently disabled. It doesn't do what the Hydrotherm version 3
! *** documentation says on p. 5-28
              WRITE (fustdout, '(A)') '*** ERROR - Option no. 22 for calculating temperature '//  &
                 'distribution is currently disabled. Please use an alternative option.  ***'
              WRITE (fuclog, '(A)') '*** ERROR - Option no. 22 for calculating temperature '//  &
                 'distribution is currently disabled. Please use an alternative option.  ***'
              ierr(173) = .TRUE.
              DEALLOCATE(dumxy)
              RETURN
!
!              kod(10) = 1
!              kodt = 2                        ! ... temperature input for subset of nodes
!              ncols = 0
!              prflag = -1                     ! ... turn off b.c. h,T printing
!              DO 
!281              READ(fuins,*) i, j, temptopalt, tempbotalt
!                 IF(i == 0 .AND. j == 0) EXIT     ! ... Terminate specified column input 
!                 ! ...            check that I,J for column is OK
!                 IF (i > nx .OR. j > ny) THEN
!                    WRITE (fustdout,9015) 'Temperature', i, j
!                    WRITE (fuclog,9015) 'Temperature', i, j
9015                FORMAT (/tr5,'**** INPUT ERROR: ', a, ' ***'/tr10,  &
                         'Column designation is invalid for I,J= ',2I4)
!                    ierr(52) = .TRUE.
!                    GO TO 281
!                 END IF
!                 ncols = ncols + 1             ! ... count the number of columns with data input
!                 IF (temptopalt > 0.) temptop = temptopalt
!                 IF (tempbotalt > 0.) tempbot = tempbotalt
!                 icup = 0
!                 DO  k = nz, 1, -1
!                    ic = (k-1)*nx + i + (j-1)*nx*nz
!                    ik = (k-1)*nx + i
!                    ! ...             skip past nodes not in the domain
!                    IF (np(ik,j) == 0) THEN
!                       tcassoc(ic) = 0._kdp
!                       prflag(ik,j) = -1          ! ... no print at this node
!                       CYCLE
!                    END IF
!                    ! ... assign input temperature to top cell in domain in each column
!                    IF (icup == 0) THEN
!                       tcassoc(ic) = temptop
!                       ehassoc(ic) = -99._kdp
!                       prflag(ik,j) = 1          ! ... print at this node
!                    ELSE
!                       ! ... calc temperature at lower nodes using input thermal gradient
!                       tcassoc(ic) = tcassoc(icup) + tempbot*0.5_kdp*(dz(k+1)+dz(k))/1.e5_kdp
!                       ehassoc(ic) = -99._kdp
!                       prflag(ik,j) = 1          ! ... print at this node
!                    END IF
!                    icup = ic
!                 END DO
!              END DO
!              IF (ncols == nxy) kodt = 1         ! ... temperature data have been input for all columns
!              ! ...        check if data entered for more columns than exist in model
!              IF (ncols > nxy) THEN
!                 WRITE (fustdout,9020) 'Temperature', ncols, nxy
!                 WRITE (fuclog,9020) 'Temperature', ncols, nxy
9020             FORMAT (/tr5,'**** INPUT ERROR: ',a,' ***'/tr10,  &
                      'Data input for more columns (',i4,') than exist in mesh (',i4,')')
!                 ierr(53) = .TRUE.
!              END IF
           ELSE IF (iopttemp == 3) THEN
              ! ... Temperature option 3:  boiling point with depth curve
              ! ...             hydrostatic pressure gradient must be requested
              ! ...             subroutine PRESBOIL calculates the boiling point curve in INPUTPH
              kodt = 3           ! ... temperature data input for boiling point options
              kod(10) = 70
           ELSE IF (iopttemp == 33) THEN
              ! ... Temperature option 33:  identical to 3 except that the
              ! ...        boiling point with depth curve is only for specified columns (I,J)
! *** This option is currently disabled. It doesn't do what the Hydrotherm version 3
! *** documentation says on p. 5-29
              WRITE (fustdout, '(A)') '*** ERROR - Option no. 33 for calculating temperature '//  &
                 'distribution is currently disabled. Please use an alternative option.  ***'
              WRITE (fuclog, '(A)') '*** ERROR - Option no. 33 for calculating temperature '//  &
                 'distribution is currently disabled. Please use an alternative option.  ***'
              ierr(174) = .TRUE.
              DEALLOCATE(dumxy)
              RETURN
!
!              kodt = 3
!              kod(10) = 75
!              DO 
!                 READ(fuins,*) i, j, ptop
!                 IF(i == 0 .AND. j == 0) EXIT     ! ... Terminate specified column input
!                 ! ...      check that I,J for column is OK
!                 IF (i > nx .OR. j > ny) THEN
!                    WRITE (fustdout,9015) 'Temperature', i, j
!                    WRITE (fuclog,9015) 'Temperature', i, j
!                    ierr(55) = .TRUE.
!                 END IF
!                 ! ...     check that PTOP is within limits of table
!                 IF (ptop > 1.e10_kdp) THEN
!                    !... special patch for low pressures
!                    !$$            IF (Ptop.GT.1.0D10 .OR. Ptop.LT.0.5D6) THEN
!                    WRITE (fustdout,9055) ptop
!                    WRITE (fuclog,9055) ptop
!9055                FORMAT (' **** INPUT ERROR: P at top of region outside of table limits ***', /tr10,  &
!                         ' P = ', 1pg12.5)
!                    ierr(54) = .TRUE.
!                 END IF
!                 CALL presboil(p, ehassoc, tcassoc, i, j)
!                 ! ... Enthalpy is also calculated here
!              END DO
           ELSE
              GO TO 990
           END IF
        ELSE IF (array == pres) THEN           ! ...                 ----- Pressure -----
           ! ... Read in Pressure at the top of the region mesh and compute
           ! ...     hydrostatic pressure distribution in subroutine INPUTPH
           IF (informt == 1) THEN
              READ(fuins,*,iostat=rd_err) ptop     ! ... old format
              ioptpres = 1
           ELSE
              READ(fuins,*,iostat=rd_err) ioptpres, ptop
           END IF
           IF(rd_err /=0) THEN
              WRITE (fustdout,'(/a)') '**** INPUT ERROR: Data for CALC option for Pressure *** '
              WRITE (fuclog,'(/a)') '**** INPUT ERROR: Data for CALC option for Pressure *** '
              ierr(60) = .TRUE.
              DEALLOCATE(dumxy)
              RETURN
           END IF
           IF (ioptpres == 1) THEN
              ! ... Pressure option 1: hydrostatic pressure gradient from top pressure
              ! ...      uniform over entire mesh
              kod(11) = 60
              DO ij=1,nxy
                 ptopa(ij) = ptop
              END DO
           ELSE IF (ioptpres == 11) THEN
              ! ... Pressure option 11: piecewise linear pressure profile
              ! ...      from top pressure
              ! ...      PTOP to bottom pressure PBOT with up to 5 linear segments
              kod(11) = 3                 ! ...  Set KOD(11)=3 for uniform row values
              IF(nx == 1) kod(11) = 1     ! ...      except when NX=1  set KOD(11)=1
              ! ... z values are input from top to bottom as depth values and are 
              ! ...      indexed to relate to k increasing from
              ! ...      bottom to top of region. Units match z-coordinate input units
              READ(fuins,*) nzppts,(zp(i),pvz(i),i=nzppts,1,-1)
              IF(nzppts > 6) THEN
                 ierr(67) = .TRUE.
              END IF
              ! ... Convert from depth to z-coordinate and into cm
              DO  i=1,nzppts
                 zp(i) = zz(nz)+0.5_kdp*dz(nz) - zp(i)*unitfac(14)
              END DO
              DO  k=1,nz
                 ! ... interpolate and convert to cgs units for pressure
                 presstemp = interp1(zz(k),nzppts,zp,pvz)*unitfac(11)
                 DO j = 1,ny
                    DO i=1,nx
                       ik = (k-1)*nx + i
                       ic = (k-1)*nx + i + (j-1)*nxz
                       ! ... Only load the active cells
                       IF (np(ik,j) == 0) THEN
                          p(ic) = 0._kdp
                       ELSE
                          p(ic) = presstemp
                       END IF
                       ! ... Load ptopa for later range check for confined flow
                       ij = (j-1)*nx + i
                       IF(k == ktop(ij)) ptopa(ij) = presstemp
                    END DO
                 END DO
              END DO
           ELSE IF (ioptpres == 21) THEN
              ! ... Pressure option 21: Similar to option 1 except that the
              ! ...      upper pressure values are not uniform. Pressure values must be
              ! ...      input for all columns in mesh. These values are assigned to the
              ! ...      uppermost active nodes.
              kod(11) = 60
              READ(fuins,*) (ptopa(ij),ij=1,nxy)
           ELSE IF (ioptpres == 22) THEN
              ! ... Pressure option 22: Similar to option 1 except that the
              ! ...      top pressure is a default value with other top pressures only
              ! ...      for specified columns (I,J).
              kod(11) = 60
              ! ... Assign default pressure
              DO ij=1,nxy
                 ptopa(ij) = ptop
              END DO
              ncols = 0
              DO
                 READ(fuins,*) i, j, ptopalt
                 IF (i == 0) EXIT
                 ! ... Check that I,J for column is valid
                 IF (i < 0 .OR. i > nx .OR. j < 0 .OR. j > ny) THEN
                    WRITE (fustdout,9015) 'pressure', i, j
                    WRITE (fuclog,9015) 'pressure', i, j
                    ierr(130) = .TRUE.
                 ELSE
                    ij = (j-1)*nx + i
                    IF (ptopalt > 0.) ptopa(ij) = ptopalt
                    ncols = ncols + 1      ! ... Count the number of columns with data input
                 END IF
              END DO       ! ... end of do forever
              ! ... Check if data have been entered for more columns than exist in mesh
              IF(ncols > nxy) THEN
                 WRITE (fustdout,9020) 'pressure', ncols, nxy
                 WRITE (fuclog,9020) 'pressure', ncols, nxy
                 ierr(131) = .TRUE.
              END IF
           END IF
           ! ... Convert from user specified units to base units (cgs)
           ptopa = ptopa*unitfac(11)
           ! ... Check that PTOPA is within limits of the enthalpy tables
           DO ij=1,nxy
              IF (ptopa(ij) > 1.0e10_kdp .OR. (.NOT.unconfined .AND. ptopa(ij) < 0.5e6_kdp)) THEN
                 CALL mtoijk(ij,i,j,k,nx,ny)
                 WRITE (fustdout, '(A/tr10,A,1PG12.5,a,2i4)')  &
                      '**** INPUT ERROR: P at top of mesh is beyond table limits ***',' P=  ',ptopa(ij),  &
                      ' for I,J= ',i,j
                 WRITE (fuclog, '(A/tr10,A,1PG12.5,a,2i4)')  &
                      '**** INPUT ERROR: P at top of mesh is beyond table limits ***',' P=  ',ptopa(ij),  &
                      ' for I,J= ',i,j
                 ierr(56) = .TRUE.
              END IF
           END DO
        ELSE IF (array == precip) THEN
           ! ... Precipitation Flux at top of the region
           ! ...     compute Precipitation flux input using a linear interpolation in the
           ! ...       X-direction.  In WELLALLO the flux (cm^3/s-cm^2) is multiplied by the
           ! ...       area (cm^2) of each cell to give a flow rate (cm^3/s)
           READ(fuins,*) icalcopt, flxleft, flxright
           ! ...         option 1 (presently the only option)
           IF (icalcopt == 1) THEN
              ! ...                if two values are equal
              IF (flxleft == flxright) THEN
                 kod(18) = 2
                 DO  ij = 1, nxy
                    qprecip(ij) = flxleft
                 END DO
              ELSE
                 ! ... check for input error using variable precipitation flux with NX=1
                 ! ... this causes division by 0
                 IF (nx == 1) THEN
                    WRITE (fustdout, '(/,2A)')  &
                         '**** INPUT ERROR: Spatially Variable Precipitation flux is not ',  &
                         'valid for NX=1 ***'
                    WRITE (fuclog, '(/,2A)')  &
                         '**** INPUT ERROR: Spatially Variable Precipitation flux is not ',  &
                         'valid for NX=1 ***'
                    ierr(57) = .TRUE.
                 END IF
                 kod(18) = 1
                 DO  i = 1, nx
                    xratio = (xx(i)-xx(1))/(xx(nx)-xx(1))
                    DO  j = 1, ny
                       ij = (j-1)*nx + i
                       qprecip(ij) = (1._kdp-xratio)*flxleft+xratio*flxright
                    END DO
                 END DO
              END IF
              ! ...            convert from user specified units to base units (cgs)
              DO  ij = 1, nxy
                 qprecip(ij) = qprecip(ij)*unitfac(18)
              END DO
              CYCLE loop310
           END IF
           GO TO 990
        ELSE IF (array == tfluxn) THEN
           ! ... Associated Temperature for Precipitation Flux at top of the region
           ! ...     compute associated temperature for Precipitation flux input using
           ! ...     linear interpolation in the
           ! ...       X-direction.  In WELLALLO the flux (cm^3/s-cm^2) is multiplied
           ! ...      by the
           ! ...       associated enthalpy and area (cm^2) of each cell to give a
           ! ...      heat flow rate (erg/s)
           READ(fuins,*) icalcopt, tflxleft, tflxright
           ! ...         option 1 (presently the only option)
           IF (icalcopt == 1) THEN
              IF (tflxleft == tflxright) THEN     ! ... if the two values are equal
                 kod(19) = 2
                 DO  ij=1,nxy
                    tflux(ij) = tflxleft
                 END DO
              ELSE
                 ! ...  input error if using variable precipitation flux with NX=1
                 IF (nx == 1) THEN
                    WRITE (fustdout, '(/2A)')  &
                         '**** INPUT ERROR: Spatially Variable associated temperature is not',  &
                         ' valid for NX=1 ***'
                    WRITE (fuclog, '(/2A)')  &
                         '**** INPUT ERROR: Spatially Variable associated temperature is not',  &
                         ' valid for NX=1 ***'
                    ierr(58) = .TRUE.
                    DEALLOCATE(dumxy)
                    RETURN
                 END IF
                 kod(19) = 1
                 DO  i=1,nx
                    xratio = (xx(i)-xx(1))/(xx(nx)-xx(1))
                    DO  j=1,ny
                       ij = (j-1)*nx + i
                       tflux(ij) = (1._kdp-xratio)*tflxleft+xratio*tflxright
                    END DO
                 END DO
              END IF
              CYCLE loop310
           END IF
           GO TO 990
        ELSE
           GO TO 990     ! ... write an error message if CALC was selected for
           ! ...               an array for which no function is yet available
        END IF
     ELSE IF (duminput(1:8) == 'RTFDEPTH') THEN
        ! ...  RTFDEPTH - by rock type with functions of depth
        IF (.NOT.ALLOCATED(zone_prop)) ALLOCATE(zone_prop(nrxtype),  &
             STAT=a_err)
        IF (a_err /= 0) THEN  
           PRINT *, "Array allocation failed: readarry, zone_prop"
           ierr(198) = .TRUE.
           DEALLOCATE(dumxy)
           RETURN
        ENDIF
        IF (array == xper) THEN                ! ...  X Permeability
           kod(2) = 1     ! ... variable value array
           READ(fuins,*) dep_sca_opt, z_datum
           DO iz=1,nrxtype
              READ(fuins,*) zone_no, zone_prop(zone_no)%fn_no
              IF(zone_prop(zone_no)%fn_no == 1) THEN
                 READ(fuins,*) zone_prop(zone_no)%par1
              ELSE 
                 READ(fuins,*) zone_prop(zone_no)%par1, zone_prop(zone_no)%par2,  &
                      zone_prop(zone_no)%par3, zone_prop(zone_no)%par4
              END IF
           END DO
           ! ... Calculate the permeability distribution value at each node
           intyp = dep_sca_opt
           z_datum = z_datum*unitfac(14)   ! ... Convert from user specified units to base units (cgs)
           DO j=1,ny
              DO k=1,nz
                 DO i=1,nx
                    ij = (j-1)*nx + i
                    ic = (k-1)*nx + i + (j-1)*nxz
                    zone_no = icrxtype(ic)
                    ifunct = zone_prop(zone_no)%fn_no
                    IF (ilrxtype(zone_no)) THEN               ! ... a legal zone number
                       IF (ifunct == 1) THEN
                          rxparm(2,zone_no) = -16._kdp
                          ! ... Permeability option 1: constant function of z-elevation or depth from 
                          ! ...      z-datum or depth from land surface
                          ! ...      uniform in horizontal directions
                          ! ... Convert from user specified units to base units (cgs)
                          perm1 = zone_prop(zone_no)%par1*unitfac(2)
                          xk(ic) = perm1
                       ELSE IF (ifunct == 2) THEN
                          rxparm(2,zone_no) = -32._kdp
                          ! ... Permeability option 2: linear function of z-elevation or depth from 
                          ! ...      z-datum or depth from land surface
                          ! ...      uniform in horizontal directions
                          IF (intyp == 1) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             zele1 = zone_prop(zone_no)%par1*unitfac(14)
                             zele2 = zone_prop(zone_no)%par3*unitfac(14)
                             perm1 = zone_prop(zone_no)%par2*unitfac(2)
                             perm2 = zone_prop(zone_no)%par4*unitfac(2)
                             xk(ic) = perm1 + (perm2-perm1)*(zz(k)-zele1)/(zele2-zele1)
                          ELSE IF (intyp == 2) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             depth1 = zone_prop(zone_no)%par1*unitfac(14)
                             depth2 = zone_prop(zone_no)%par3*unitfac(14)
                             perm1 = zone_prop(zone_no)%par2*unitfac(2)
                             perm2 = zone_prop(zone_no)%par4*unitfac(2)
                             xk(ic) = perm1 + (perm2-perm1)*(z_datum-zz(k)-depth1)/(depth2-depth1)
                          ELSE IF (intyp == 3) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             depth1 = zone_prop(zone_no)%par1*unitfac(14)
                             depth2 = zone_prop(zone_no)%par3*unitfac(14)
                             perm1 = zone_prop(zone_no)%par2*unitfac(2)
                             perm2 = zone_prop(zone_no)%par4*unitfac(2)
                             xk(ic) = perm1 + (perm2-perm1)*(zls(ij)-zz(k)-depth1)/(depth2-depth1)
                          END IF
                       ELSE IF (ifunct == 3) THEN
                          rxparm(2,zone_no) = -256._kdp
                          ! ... Permeability option 3: linear variation of logarithm of permeability      
                          ! ...      with z-elevation or depth from z-datum or depth from 
                          ! ...      land surface, using log base 10
                          ! ...      from Williams and Narasimhan (1989)
                          ! ...      uniform in horizontal directions
                          IF (intyp == 1) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             zele1 = zone_prop(zone_no)%par1*unitfac(14)
                             zele2 = zone_prop(zone_no)%par3*unitfac(14)
                             perm1 = zone_prop(zone_no)%par2*unitfac(2)
                             perm2 = zone_prop(zone_no)%par4*unitfac(2)
                             a0 = (zele1*LOG10(perm2)-zele2*LOG10(perm1))/(zele1-zele2)
                             a1 = LOG10(perm1/perm2)/(zele1-zele2)
                             xk(ic) = 10._kdp**(a0+a1*zz(k))
                          ELSE IF (intyp == 2) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             depth1 = zone_prop(zone_no)%par1*unitfac(14)
                             depth2 = zone_prop(zone_no)%par3*unitfac(14)
                             zele1 = z_datum-depth1
                             zele2 = z_datum-depth2
                             perm1 = zone_prop(zone_no)%par2*unitfac(2)
                             perm2 = zone_prop(zone_no)%par4*unitfac(2)
                             a0 = (zele1*LOG10(perm2)-zele2*LOG10(perm1))/(zele1-zele2)
                             a1 = LOG10(perm1/perm2)/(zele1-zele2)
                             xk(ic) = 10**(a0+a1*zz(k))
                          ELSE IF (intyp == 3) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             depth1 = zone_prop(zone_no)%par1*unitfac(14)
                             depth2 = zone_prop(zone_no)%par3*unitfac(14)
                             perm1 = zone_prop(zone_no)%par2*unitfac(2)
                             perm2 = zone_prop(zone_no)%par4*unitfac(2)
                             a0 = (depth1*LOG10(perm2)-depth2*LOG10(perm1))/(depth1-depth2)
                             a1 = LOG10(perm1/perm2)/(depth1-depth2)
                             xk(ic) = 10._kdp**(a0+a1*(zls(ij)-zz(k)))
                          END IF
                       ELSE IF (ifunct == 4) THEN
                          rxparm(2,zone_no) = -128._kdp
                          ! ... Permeability option 4: variation of logarithm of permeability
                          ! ...      with logarithm of z-elevation or depth from z-datum or depth from 
                          ! ...      land surface, using log base 10
                          ! ...      from Manning and Ingebritsen (1999)
                          ! ...      uniform in horizontal directions
                          IF (intyp == 1) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             zele1 = zone_prop(zone_no)%par1*unitfac(14)
                             zele2 = zone_prop(zone_no)%par3*unitfac(14)
                             depth1 = z_datum-zele1
                             depth2 = z_datum-zele2
                             perm1 = zone_prop(zone_no)%par2*unitfac(2)
                             perm2 = zone_prop(zone_no)%par4*unitfac(2)
                             a0 = (LOG10(depth1)*LOG10(perm2)-LOG10(depth2)*LOG10(perm1))/  &
                                  LOG10(depth1/depth2)
                             a1 = LOG10(perm1/perm2)/LOG10(depth1/depth2)
                             xk(ic) = 10._kdp**(a0+a1*LOG10(z_datum-zz(k)))
                          ELSE IF (intyp == 2) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             depth1 = zone_prop(zone_no)%par1*unitfac(14)
                             depth2 = zone_prop(zone_no)%par3*unitfac(14)
                             perm1 = zone_prop(zone_no)%par2*unitfac(2)
                             perm2 = zone_prop(zone_no)%par4*unitfac(2)
                             a0 = (LOG10(depth1)*LOG10(perm2)-LOG10(depth2)*LOG10(perm1))/  &
                                  LOG10(depth1/depth2)
                             a1 = LOG10(perm1/perm2)/LOG10(depth1/depth2)
                             xk(ic) = 10._kdp**(a0+a1*LOG10(z_datum-zz(k)))
                          ELSE IF (intyp == 3) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             depth1 = zone_prop(zone_no)%par1*unitfac(14)
                             depth2 = zone_prop(zone_no)%par3*unitfac(14)
                             perm1 = zone_prop(zone_no)%par2*unitfac(2)
                             perm2 = zone_prop(zone_no)%par4*unitfac(2)
                             a0 = (LOG10(depth1)*LOG10(perm2)-LOG10(depth2)*LOG10(perm1))/  &
                                  LOG10(depth1/depth2)
                             a1 = LOG10(perm1/perm2)/LOG10(depth1/depth2)
                             xk(ic) = 10._kdp**(a0+a1*LOG10(zls(ij)-zz(k)))
                          END IF
                       ELSE
                          GO TO 990    ! ... invalid option
                       END IF
                    ELSE
                       ! ... invalid rock type or zone index - write error statement
                       WRITE (fustdout,9010) 'permeability'
                       WRITE (fustdout,9020) fmt(1:7)
                       WRITE (fustdout,9130) irx
                       WRITE (fuclog,9010) 'permeability'
                       WRITE (fuclog,9020) fmt(1:7)
                       WRITE (fuclog,9130) irx
9130                   FORMAT (tr10,'Rock type ', i2, ' not used in this simulation')
                       ierr(24) = .TRUE.
                    END IF
                 END DO
              END DO
           END DO
        ELSE IF (array == 'PORO') THEN               ! ... Porosity
           kod(1) = 1     ! ... variable value array
           ! ... Read the data
           READ(fuins,*) dep_sca_opt, z_datum
           DO iz=1,nrxtype
              READ(fuins,*) zone_no, zone_prop(zone_no)%fn_no
              IF(zone_prop(zone_no)%fn_no == 1) THEN
                 READ(fuins,*) zone_prop(zone_no)%par1
              ELSE 
                 READ(fuins,*) zone_prop(zone_no)%par1, zone_prop(zone_no)%par2,  &
                      zone_prop(zone_no)%par3, zone_prop(zone_no)%par4
              END IF
           END DO
           ! ... Calculate the porosity distribution value at each node
           intyp = dep_sca_opt
           z_datum = z_datum*unitfac(14)   ! ... Convert from user specified units to base units (cgs)
           DO j=1,ny
              DO k=1,nz
                 DO i=1,nx
                    ij = (j-1)*nx + i
                    ic = (k-1)*nx + i + (j-1)*nxz
                    zone_no = icrxtype(ic)
                    ifunct = zone_prop(zone_no)%fn_no
                    IF (ilrxtype(zone_no)) THEN          ! ... a legal zone number
                       IF (ifunct == 1) THEN
                          rxparm(1,zone_no) = -16._kdp
                          ! ... Porosity option 1: constant function of z-elevation or depth from 
                          ! ...      z-datum or depth from land surface
                          ! ...      uniform in horizontal directions
                          ! ... Convert from user specified units to base units (cgs)
                          poro1 = zone_prop(zone_no)%par1*unitfac(2)
                          phi(ic) = poro1
                       ELSE IF (ifunct == 2) THEN
                          rxparm(1,zone_no) = -32._kdp
                          ! ... Porosity option 2: linear function of z-elevation or depth from 
                          ! ...      z-datum or depth from land surface
                          ! ...      uniform in horizontal directions
                          IF (intyp == 1) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             zele1 = zone_prop(zone_no)%par1*unitfac(14)
                             zele2 = zone_prop(zone_no)%par3*unitfac(14)
                             poro1 = zone_prop(zone_no)%par2
                             poro2 = zone_prop(zone_no)%par4
                             phi(ic) = poro1 + (poro2-poro1)*(zz(k)-zele1)/(zele2-zele1)
                          ELSE IF (intyp == 2) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             depth1 = zone_prop(zone_no)%par1*unitfac(14)
                             depth2 = zone_prop(zone_no)%par3*unitfac(14)
                             poro1 = zone_prop(zone_no)%par2
                             poro2 = zone_prop(zone_no)%par4
                             phi(ic) = poro1 + (poro2-poro1)*(z_datum-zz(k)-depth1)/(depth2-depth1)
                          ELSE IF (intyp == 3) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             depth1 = zone_prop(zone_no)%par1*unitfac(14)
                             depth2 = zone_prop(zone_no)%par3*unitfac(14)
                             poro1 = zone_prop(zone_no)%par2
                             poro2 = zone_prop(zone_no)%par4
                             phi(ic) = poro1 + (poro2-poro1)*(zls(ij)-zz(k)-depth1)/(depth2-depth1)
                          END IF
                       ELSE IF (ifunct == 3) THEN
                          rxparm(1,zone_no) = -64._kdp
                          ! ... Porosity option 3: exponential decrease of porosity 
                          ! ...      with z-elevation or depth from z-datum or depth from 
                          ! ...      land surface, based on power of 10
                          ! ...      Equivalent to log porosity is a linear function of
                          ! ...      z-elevation or depth
                          ! ...      uniform in horizontal directions
                          IF (intyp == 1) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             zele1 = zone_prop(zone_no)%par1*unitfac(14)
                             zele2 = zone_prop(zone_no)%par3*unitfac(14)
                             poro1 = zone_prop(zone_no)%par2
                             poro2 = zone_prop(zone_no)%par4
                             a0 = (zele1*LOG10(poro2)-zele2*LOG10(poro1))/(zele1-zele2)
                             a1 = LOG10(poro1/poro2)/(zele1-zele2)
                             phi(ic) = a0*10._kdp**(a1*zz(k))
                          ELSE IF (intyp == 2) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             depth1 = zone_prop(zone_no)%par1*unitfac(14)
                             depth2 = zone_prop(zone_no)%par3*unitfac(14)
                             zele1 = z_datum-depth1
                             zele2 = z_datum-depth2
                             poro1 = zone_prop(zone_no)%par2
                             poro2 = zone_prop(zone_no)%par4
                             a0 = (zele1*LOG10(poro2)-zele2*LOG10(poro1))/(zele1-zele2)
                             a1 = LOG10(poro1/poro2)/(zele1-zele2)
                             phi(ic) = a0*10._kdp**(a1*zz(k))
                          ELSE IF (intyp == 3) THEN
                             ! ... Convert from user specified units to base units (cgs)
                             depth1 = zone_prop(zone_no)%par1*unitfac(14)
                             depth2 = zone_prop(zone_no)%par3*unitfac(14)
                             poro1 = zone_prop(zone_no)%par2
                             poro2 = zone_prop(zone_no)%par4
                             a0 = 10**(depth1*LOG10(poro2)-depth2*LOG10(poro1))/(depth1-depth2)
                             a1 = LOG10(poro1/poro2)/(depth1-depth2)
                             phi(ic) = a0*10._kdp**(a1*(zls(ij)-zz(k)))
                          END IF
                       ELSE
                          GO TO 990     ! ... invalid option
                       END IF
                    ELSE
                       ! ... invalid rock type or zone index - write error statement
                       WRITE (fustdout,9010) 'porosity'
                       WRITE (fustdout,9020) fmt(1:7)
                       WRITE (fustdout,9130) irx
                       WRITE (fuclog,9010) 'porosity'
                       WRITE (fuclog,9020) fmt(1:7)
                       WRITE (fuclog,9130) irx
                       ierr(24) = .TRUE.
                    END IF
                 END DO
              END DO
           END DO
        ELSE
           GO TO 990     ! ... write an error message if RTFDEPTH was selected for
           ! ...               an array for which no function is yet available
        END IF
        DEALLOCATE(zone_prop,  &
             STAT=da_err)
        IF (da_err /= 0) THEN  
           PRINT *, "Array deallocation failed: readarry, zone_prop"
!!$           ierr(198) = .TRUE.
           DEALLOCATE(dumxy)
           RETURN
        ENDIF
     END IF
  END DO loop310
! ... Set flag for calculating functions of temperature
  lrxftreq = .false.
  do kpar=1,mprxparm
     if(irxftopt(kpar,1,1) > 1) lrxftreq = .true.
  end do
  lrxftimreq = .false.
  do kpar=2,4
     if(irxftimopt(kpar,1,1) > 1) lrxftimreq = .true.
  end do
  DEALLOCATE(dumxy)
  RETURN
  ! ############## write error message; improper format selected ##############
990 WRITE (fustdout,9010) fmt(1:7), array
  WRITE (fuclog,9010) fmt(1:7), array
9010 FORMAT ('**** INPUT ERROR: ', a7, ' format not valid for ', a4, ' ***')
  ierr(59) = .TRUE.
  DEALLOCATE(dumxy)
  RETURN
995 WRITE (fustdout,9035)  '**** INPUT ERROR: Problem with RXKX input *****',  &
       'Check that', nrxtype,' rock types are defined'
  WRITE (fuclog,9035) '**** INPUT ERROR: Problem with RXKX input *****',  &
       'Check that', nrxtype,' rock types are defined'
9035 FORMAT (/tr5,a/tr9,a,i2,a)
  ierr(33) = .TRUE.
  DEALLOCATE(dumxy)
END SUBROUTINE readarry
