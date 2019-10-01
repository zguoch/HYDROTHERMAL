SUBROUTINE init_wr_pmp(initial,flag)
  ! ... Purpose: To set up and write the initial or new values for the porous
  ! ...      matrix properties - porosity, permeability, thermal conductivity,
  ! ...      specific heat, rock density, compressibility.
  ! ... Includes functions of rock type, temperature, and time
  ! ... Functions of depth do not have their parameters stored. Instead, the
  ! ...     initial spatial distributions of the relevant properties are retained.
  USE f_units
  USE bc, ONLY: ibc
  USE control
  USE i_c
  USE mesh
  USE parameters
  USE units
  USE variables
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: initial       ! ... TRUE - initial values input
                                       ! ... FALSE - new values input
  INTEGER, INTENT(INOUT) :: flag
  !
  CHARACTER (LEN=80) :: ufile
  CHARACTER (LEN=21), DIMENSION(mprxparm) :: chrrxparm=(/  &
       'Porosity:            ', 'X-Permeability:      ',  &
       'Y-Permeability:      ','Z-Permeability:      ', 'Thermal Conductivity:',  &
       'Specific Heat:       ','Rock Density:        ', 'Compressibility:     '/)
  CHARACTER (LEN=10), DIMENSION(6) :: chstep=(/  &
       'Parameter ',' for      ','   T <    ', ' <= T <   ',' <= T     ','Factor    '/)
  CHARACTER (LEN=24), DIMENSION(3) :: chline=(/'Parameter, Temperature  ',  &
       '                        ','Factor, Temperature     '/)
  CHARACTER (LEN=16), DIMENSION(3) :: chline2=(/'Parameter, Time ','                ',  &
       'Factor, Time    '/)
  CHARACTER(LEN=80) :: cfmta, cfmtb, cfmtc
  CHARACTER(LEN=42), PARAMETER :: cfmta1="(/'Column (I)->',i3,13I5/(12X,13I5))    ",  &
      cfmta2="(/'Column (I)->',i3,23I5/(12X,23I5))   ",  &
      cfmta3="(/'Column (I)->',i3,28I5/(12X,28I5))   ",  &
      cfmtb1="('Row',i3,' -> ',14I5/(12X,13I5))      ",  &
      cfmtb2="('Row',i3,' -> ',24I5/(12X,23I5))      ",  &
      cfmtb3="('Row',i3,' -> ',29I5/(12X,28I5))      ",  &
      cfmtc1="(1P,10(g11.5,1X)/(10(g11.5,1X)))       ",  &
      cfmtc2="(1P,10(g11.5,1X))                      "
  INTEGER :: i, ic, icxf, icxl, irx, j, k, kkk, kodd, kodparm, lll
  INTEGER :: funct_no
  REAL(KIND=kdp) :: time_pr
  INTERFACE
     SUBROUTINE phreg(p,h)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
     END SUBROUTINE phreg
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
     SUBROUTINE wrnew(funit,dumarray,kodd,convfac)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: funit
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: dumarray
       INTEGER, INTENT(IN) :: kodd
       REAL(KIND=kdp), INTENT(IN OUT) :: convfac
     END SUBROUTINE wrnew
     SUBROUTINE tempdegc(p,h,tc)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
       REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: tc
     END SUBROUTINE tempdegc
  END INTERFACE
  !     ------------------------------------------------------------------
  ! ...   Kod(kodparm) =1 - array input as variable value array
  ! ...             =2 - array input as constant value array
  ! ...             =3 - array input as values constant for each row
  ! ...             =4 - array input as values constant for each column
  ! ...             5 - array input as rocktype values
  ! ...             6 - array input as rocktype with function of temperature
  ! ...             7 - array input as rocktype with function of time
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
  ! ...          1 if temperature data are input for all nodes
  ! ...          2 if temperature data are input for only a few nodes
  ! ...          3 if temperature data are input with boiling point options
  !
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.8 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  ! ... calculate permeability as f(Temperature) if required
  IF (lrxftreq) THEN
     ! ...        determine phase region
     CALL phreg(p,h)
     CALL tempdegc(p,h,tc)     ! ... calculate Temperature for permeability calculations
     CALL rockparm
  END IF
  IF (lrxftimreq) THEN
     CALL permftime            ! ... calculate permeability as f(time) if required
  END IF
  ! ...          ----- Write Rock Type Distribution and parameters ----
  ! ...                if rock type format is used
  IF (ifmt == 2) THEN
     IF (ioptpr(2) == 0) THEN
        cfmta = cfmta1
        cfmtb = cfmtb1
     ELSE IF (ioptpr(2) == 1) THEN
        cfmta = cfmta2
        cfmtb = cfmtb2
     ELSE
        cfmta = cfmta3
        cfmtb = cfmtb3
     END IF
     WRITE(fupdef,'(/A)') '--- Rock-type Distribution ---'
     WRITE(fupdef,'(10X,A,A)') '     ', 'positive integer = Active node;  0 = Inactive node, out of domain'
     DO  j = 1,ny
        IF (ny > 1) WRITE(fupdef,'(/tr7,a,i4)') 'Rock-type Numbers for Slice', j
        WRITE(fupdef, cfmta) (i,i=1,nx)
        WRITE(fupdef, '(A)') ' '
        DO  k = nz,1,-1
           icxf = (k-1)*nx + 1 + (j-1)*nx*nz
           icxl = icxf + nx - 1
           WRITE(fupdef,cfmtb) k,(icrxtype(ic), ic=icxf,icxl)
        END DO
     END DO
     WRITE(fupdef, '(/A)') '--- Parameter Values for Rock Types ---'
     DO  lll = 1,nrxtype
        irx = irxused(lll)
        WRITE(fupdef,'(/tr7,a,i3)') 'Rock Type ', irx
        DO  kodparm = 1,mprxparm
           IF (rxparm(kodparm,irx) >= 0.0_kdp) THEN
              WRITE(fupdef,9005) chrrxparm(kodparm),  &
                   rxparm(kodparm,irx)/unitfac(kodparm), unitlabel(kodparm)
9005          FORMAT (tr10, a,'           ',1PG10.3,tr1,a)
           ELSE IF (rxparm(kodparm,irx) == -99.999_kdp) THEN
              WRITE(fupdef, 9110) chrrxparm(kodparm), unitlabel(kodparm),  &
                   '- Value not defined by rock type'
9110          FORMAT (tr10,3a)
           ELSE IF (rxparm(kodparm,irx) == -16._kdp) THEN
              WRITE(fupdef, 9110) chrrxparm(kodparm), unitlabel(kodparm),  &
                   '- Value defined as uniform f(depth)'
           ELSE IF (rxparm(kodparm,irx) == -32._kdp) THEN
              WRITE(fupdef, 9110) chrrxparm(kodparm), unitlabel(kodparm),  &
                   '- Value defined as linear f(depth)'
           ELSE IF (rxparm(kodparm,irx) == -64._kdp) THEN
              WRITE(fupdef, 9110) chrrxparm(kodparm), unitlabel(kodparm),  &
                   '- Value defined as exponential f(depth)'
           ELSE IF (rxparm(kodparm,irx) == -128._kdp) THEN
              WRITE(fupdef, 9110) chrrxparm(kodparm), unitlabel(kodparm),  &
                   '- Log Value defined as f(log(depth below land surface))'
           ELSE IF (rxparm(kodparm,irx) == -256._kdp) THEN
              WRITE(fupdef, 9110) chrrxparm(kodparm), unitlabel(kodparm),  &
                   '- Log Value defined as linear f(depth below land surface)'
           ELSE
              funct_no = irxftopt(kodparm,irx,1)
              ! ... RTFTemperature parameter - option 2  piecewise constant (i.e. step function)
              IF (funct_no == 2) THEN
                 WRITE(fupdef, 9015)  &
                      chrrxparm(kodparm), unitlabel(kodparm),  &
                      chstep(1), rxftparm(kodparm,irx,1)/unitfac(kodparm),  &
                      chstep(2), chstep(3), rxftparm(kodparm,irx,2),  &
                      (chstep(1), rxftparm(kodparm,irx,2*kkk-1)/unitfac(kodparm),  &
                      chstep(2), rxftparm(kodparm,irx,2*(kkk-1)),  &
                      chstep(4), rxftparm(kodparm,irx,2*kkk),  &
                      kkk=2,irxftopt(kodparm,irx,2)-1),  &
                      chstep(1), rxftparm(kodparm,irx,2*irxftopt(kodparm,irx,2)-1)/unitfac(kodparm),  &
                      chstep(2), rxftparm(kodparm,irx,2*(irxftopt(kodparm,irx,2)-1)),  &
                      chstep(5)
9015             FORMAT (tr10, a, a, '- Value defined by Step f(Temperature)'/  &
                       tr15,a10,1pg10.3,a5,tr7,a8,0pf7.2/  &
                      (tr15,a10,1pg10.3,a,0pf7.2,a8,0pf7.2))
              ELSE IF (funct_no == 3) THEN
                 WRITE(fupdef, 9020)  &
                      chrrxparm(kodparm), unitlabel(kodparm),  &
                      chline(1), rxftparm(kodparm,irx,1)/unitfac(kodparm), &
                      rxftparm(kodparm,irx,2), &
                      (chline(2), rxftparm(kodparm,irx,2*kkk-1)/unitfac(kodparm),  &
                      rxftparm(kodparm,irx,2*kkk), kkk=2,irxftopt(kodparm,irx,2))
9020             FORMAT (tr10, a, a, '- Value defined by Log linear f(Temperature)'/  &
                      (tr15,a,1Pg10.3,',  ',0Pf10.2))
              ELSE IF (funct_no == 4) THEN
                 WRITE(fupdef, 9025)  &
                      chrrxparm(kodparm), unitlabel(kodparm),  &
                      chline(1), rxftparm(kodparm,irx,1)/unitfac(kodparm),  &
                      rxftparm(kodparm,irx,2),  &
                      (chline(2), rxftparm(kodparm,irx,2*kkk-1)/unitfac(kodparm),  &
                      rxftparm(kodparm,irx,2*kkk), kkk=2,irxftopt(kodparm,irx,2))
9025             FORMAT (tr10, a, a, '- Value defined by Linear f(Temperature)'/  &
                      (tr15,a,1Pg10.3,',  ',0Pf10.2))
              ! ... Y-permeablity when X-perm is f(Temperature)
                 IF (kodparm == 3 .AND. irxftopt(2,irx,1) >= 2)  &
                     WRITE(fupdef, 9030) chrrxparm(kodparm), unitlabel(kodparm), ykrxkxfac(irx)
              ! ... Z-permeablity when X-perm is f(Temperature)
                 IF (kodparm == 4 .AND. irxftopt(2,irx,1) >= 2)  &
                     WRITE(fupdef, 9030) chrrxparm(kodparm), unitlabel(kodparm), zkrxkxfac(irx)
9030             FORMAT (10X, a, a, '- Vary from X-Perm by factor =', 1pg9.2)
              END IF
           END IF
!
           IF (kodparm >= 2 .AND. kodparm <=4) THEN
              ! ... RTFTIME parameter - option 2  linear f(Time)
              IF (irxftimopt(kodparm,irx,1) == 2) THEN
                 IF(iyr == 2) THEN
                    time_pr = year(rxftimparm(kodparm,irx,2*kkk))
                 ELSE
                    time_pr = rxftimparm(kodparm,irx,2*kkk)
                 END IF
                 WRITE(fupdef,9125)  &
                      chrrxparm(kodparm), unitlabel(kodparm),  &
                      (chline2(kkk), rxftimparm(kodparm,irx,2*kkk-1)/unitfac(kodparm),   &
                      time_pr, kkk=1, irxftimopt(kodparm,irx,2))                         
9125             FORMAT (tr10,a,a,'- Value defined by Linear f(Time)'/  &
                      (tr15,a,tr2,1Pg10.3,',',1Pg12.4))
              END IF
              ! ... Y-permeablity when X-perm is f(Time)
              IF (kodparm == 3 .AND. irxftimopt(2,irx,1) >= 2)  &
                   WRITE(fupdef,9130) chrrxparm(kodparm), unitlabel(kodparm),  &
                   ykrxkxfac(irx)
              ! ... Z-permeablity when X-perm is f(Time)
              IF (kodparm == 4 .AND. irxftimopt(2,irx,1) >= 2)  &
                   WRITE(fupdef,9130) chrrxparm(kodparm), unitlabel(kodparm),  &
                   zkrxkxfac(irx)
9130          FORMAT (tr10,a,a,'- Vary from X-Perm by factor =',1Pg9.2)
           END IF
        END DO
     END DO
  END IF
  ! ...                       POROSITY
  IF (kod(1) > 0) THEN
     kodd = 1
     IF (initial) THEN
        IF (initphi) THEN
           WRITE(fupdef, 9035) 'Initial Porosity', unitlabel(kodd)
           WRITE(fupdef, '(5X,A,A)') 'Values read from file ',TRIM(phifile)
           IF(kod(1) == 1 .OR. kod(1) == 5) THEN
              CALL printar(2,phiinit,ibc,fupdef,unitfac(1),14,000)
           ELSE
              CALL wrnew(fupdef, phiinit, kod(1), unitfac(1))
           END IF
           WRITE(fupdef,9035) 'Porosity at start of run', unitlabel(kodd)
           IF(kod(1) == 1 .OR. kod(1) == 5) THEN
              CALL printar(2,phi,ibc,fupdef,unitfac(1),14,000)
           ELSE
              CALL wrnew(fupdef, phi, kod(1), unitfac(1))
           END IF
        ELSE
           WRITE(fupdef, 9035) 'Initial Porosity', unitlabel(kodd)
           IF(kod(1) == 1 .OR. kod(1) == 5) THEN
              CALL printar(2,phi,ibc,fupdef,unitfac(1),14,000)
           ELSE
              CALL wrnew(fupdef, phi, kod(1), unitfac(1))
           END IF
        END IF
     ELSE
        WRITE(fupdef, 9035) 'New values for Porosity', unitlabel(kodd)
           IF(kod(1) == 1 .OR. kod(1) == 5) THEN
              CALL printar(2,phi,ibc,fupdef,unitfac(1),14,000)
           ELSE
              CALL wrnew(fupdef, phi, kod(1), unitfac(1))
           END IF
     END IF
  END IF
  ! ...       save the new porosity and pressure as the initial values
  ! ...       if the rock compressibility (an input value) is > 0
  ! ...          then write the new initial values to a file which can
  ! ...          be later used to restart a stopped run
  IF (.NOT.initphi) THEN
     DO  ic = 1, nxyz
        pinit(ic) = p(ic)
        phiinit(ic) = phi(ic)
     END DO
     IF (beta(1) /= 0._kdp) THEN
        ufile = ' '
        CALL gfiles(4, ufile, flag)
        IF(flag > 0) RETURN         ! ... error abort
        IF (nx > 10) THEN
           cfmtc = cfmtc1
        ELSE
           cfmtc = cfmtc2
        END IF
        WRITE(fuicpp,'(a)') 'INITIAL PRESSURE'
        DO  j = 1, ny
           DO  k = 1, nz
              icxf = (k-1)*nx + 1 + (j-1)*nx*nz
              icxl = icxf + nx - 1
              WRITE(fuicpp,cfmtc) (pinit(ic), ic=icxf, icxl)
           END DO
        END DO
        WRITE(fuicpp,'(a)') 'INITIAL POROSITY'
        DO  j = 1, ny
           DO  k = 1, nz
              icxf = (k-1)*nx + 1 + (j-1)*nx*nz
              icxl = icxf + nx - 1
              WRITE(fuicpp,cfmtc) (phiinit(ic), ic=icxf, icxl)
           END DO
        END DO
     END IF
  END IF
  ! ...                    X-permeabilities
  IF (kod(2) > 0) THEN
     kodd = 2
     IF (initial) THEN
        WRITE(fupdef, 9035) 'X-Permeability', unitlabel(kodd)
     ELSE
        WRITE(fupdef, 9035) 'New values for X-Permeability', unitlabel(kodd)
     END IF
     IF(kod(2) == 1 .OR. kod(2) == 5) THEN
        CALL printar(2,xk,ibc,fupdef,unitfac(kodd),24,000)
     ELSE
        CALL wrnew(fupdef, xk, kod(kodd), unitfac(kodd))
     END IF
  END IF
  ! ...                    Y-permeabilities
  IF (kod(3) > 0) THEN
     kodd = 3
     IF (initial) THEN
        WRITE(fupdef, 9035) 'Y-Permeability', unitlabel(kodd)
     ELSE
        WRITE(fupdef, 9035) 'New values for Y-Permeability', unitlabel(kodd)
     END IF
     IF (kod(3) /= 90) THEN
        IF(kod(3) == 1 .OR. kod(3) == 5) THEN
           CALL printar(2,yk,ibc,fupdef,unitfac(kodd),24,000)
        ELSE
           CALL wrnew(fupdef, yk, kod(kodd), unitfac(kodd))
        END IF
     ELSE IF (ykfactor == 1.0_kdp) THEN
        WRITE(fupdef,'(/10X,A)') 'Values are the same as X-Permeability'
     ELSE
        WRITE(fupdef, '(/TR10,A,1PG9.3,A)')  &
             'Values are ', ykfactor,' times X-Permeability'
     END IF
  END IF
  ! ...                    Z-permeabilities
  IF (kod(4) > 0) THEN
     kodd = 4
     IF (initial) THEN
        WRITE(fupdef, 9035) 'Z-Permeability', unitlabel(kodd)
     ELSE
        WRITE(fupdef, 9035) 'New values for Z-Permeability', unitlabel(kodd)
     END IF
     IF (kod(4) /= 90) THEN
        IF(kod(4) == 1 .OR. kod(4) == 5) THEN
           CALL printar(2,zk,ibc,fupdef,unitfac(kodd),24,000)
        ELSE
           CALL wrnew(fupdef, zk, kod(kodd), unitfac(kodd))
        END IF
     ELSE IF (zkfactor == 1.0_kdp) THEN
        WRITE(fupdef,'(/10X,A)') 'Values are the same as X-Permeability'
     ELSE
        WRITE(fupdef, '(/TR10,A,1PG9.3,A)')  &
             'Values are ', zkfactor,' times X-Permeability'
     END IF
  END IF
  ! ...                   MEDIUM THERMAL CONDUCTIVITY
  IF (kod(5) > 0) THEN
     kodd = 5
     IF (initial) THEN
        WRITE(fupdef, 9035) 'Thermal Conductivity', unitlabel(kodd)
     ELSE
        WRITE(fupdef, 9035) 'New values for Thermal Conductivity', unitlabel(kodd)
     END IF
     IF(kod(5) == 1 .OR. kod(5) == 5) THEN
        CALL printar(2,xkc,ibc,fupdef,unitfac(kodd),24,000)
     ELSE
        CALL wrnew(fupdef, xkc, kod(kodd), unitfac(kodd))
     END IF
  END IF
  ! ...                   SPECIFIC HEAT OF ROCK
  IF (kod(6) > 0) THEN
     kodd = 6
     IF (initial) THEN
        WRITE(fupdef, 9035) 'Specific Heat of Rock', unitlabel(kodd)
     ELSE
        WRITE(fupdef, 9035) 'New values for Specific Heat of Rock',  &
             unitlabel(kodd)
     END IF
     IF(kod(6) == 1 .OR. kod(6) == 5) THEN
        CALL printar(2,phfwt,ibc,fupdef,unitfac(kodd),24,000)
     ELSE
        CALL wrnew(fupdef,phfwt,kod(kodd),unitfac(kodd))
     END IF
  END IF
  ! ...                   ROCK DENSITY
  IF (kod(7) > 0) THEN
     kodd = 7
     IF (initial) THEN
        WRITE(fupdef, 9035) 'Rock Density', unitlabel(kodd)
     ELSE
        WRITE(fupdef, 9035) 'New values for Rock Density', unitlabel(kodd)
     END IF
     IF(kod(7) == 1 .OR. kod(7) == 5) THEN
        CALL printar(2,df,ibc,fupdef,unitfac(kodd),24,000)
     ELSE
        CALL wrnew(fupdef, df, kod(kodd), unitfac(kodd))
     END IF
  END IF
  ! ...                   ROCK COMPRESSIBILITY
  IF (kod(8) > 0) THEN
     kodd = 8
     IF (initial) THEN
        WRITE(fupdef, 9035) 'Rock Compressibility', unitlabel(kodd)
     ELSE
        WRITE(fupdef, 9035) 'New values for Rock Compressibility', unitlabel(kodd)
     END IF
     IF(kod(8) == 1 .OR. kod(8) == 5) THEN
        CALL printar(2,beta,ibc,fupdef,unitfac(kodd),24,000)
     ELSE
        CALL wrnew(fupdef, beta, kod(kodd), unitfac(kodd))
     END IF
  END IF
9035 FORMAT (/, '--- ', a, ' ---  ', a)
END SUBROUTINE init_wr_pmp
