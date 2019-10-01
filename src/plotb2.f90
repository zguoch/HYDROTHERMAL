SUBROUTINE plotb2
  !     Purpose:  To write results for the current time-level to file for
  !           later plotting using the Basin2 plotting program B2plot.
  !           File contains X and Z fluid velocity vectors, along with
  !           pressure, enthalpy, temperature, etc. data for each node.
  !           Written for the new "sep88" format for B2plot
  !     note: IOPTPR(1) = 8 for ASCII plotfile;  = 9 for binary plotfile
  !     note: the first three field variables must be (1) x-vel (2) z-vel
  !           (3) v-zm and (4) depth bsl.
  !           ibdy signifies that boundary values (0) are not or (1) are
  !           significant.
  !     Method:  There are a number of differences between the way
  !           HYDROTHERM handles the finite-difference grid and the way that
  !           B2plot expects the input to be formatted.  These differences
  !           include: 1) the origin for the I,J,K indicies is at the
  !           bottom of the model in HYDROTHERM and at the top in Basin2,
  !           2) HYDROTHERM is a 3-d program and Basin2 is 2-d, 3) the top
  !           row of nodes in Basin2 cannot include any nodes that are
  !           not in the domain of the problem, and they are all constant
  !           pressure and temperature nodes.
  USE machine_constants, ONLY: kdp
  USE f_units, ONLY: fupltsca
  USE parameters          !!$*** for flint
  USE control
  USE mesh
!!$  USE parameters
  USE variables
  USE units
  IMPLICIT NONE
  INTEGER, PARAMETER :: nfld=12, nld=204
  CHARACTER (LEN=40), DIMENSION(nfld), SAVE :: ufield=(/ &
       'X-velocity of water (cm/yr)            ', &
       'Z-velocity of water (cm/yr)            ', &
       'Z-velocity of steam (cm/yr)            ', &
       'Depth (km)                             ', &
       'Hydraulic potential (atm)              ', &
       'Temperature (Deg.C)                    ', &
       'Pressure (atm)                         ', &
       'Porosity (-)                           ', &
       'Log x-permeability (darcy)             ', &
       'Log z-permeability (darcy)             ', &
       'Water saturation (fraction)            ', &
       'Thermal conductivity (cal/s-cm-Deg.C)  '/)
  CHARACTER (LEN=70), DIMENSION(nld) :: ubuff
  INTEGER :: i, ic, ij, ik, ikzp, j, k, m, n, nfield=12, npoly
  !       nfld  -- max. number of fields to appear in the output
  !       nfld1 -- index of last field, exclusive of those set for rock
  !                fractions
  !       nld   -- dimension of arrays xline, zline; number of points
  !                in polygon bounding the basin
  INTEGER, DIMENSION(nfld) :: ibdy=(/1,1,1,1,1,1,1,1,0,0,1,0/)
  INTEGER, DIMENSION(6) :: ifield
  LOGICAL :: lprint
  REAL(KIND=kdp) :: width, xleft, xposition, xright, ysl, zposition, ztop
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: polyx, polyz
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: field
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: satnw2
  INTEGER :: alloc_stat
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.8 $//$Date: 2007/12/20 20:32:47 $'
  ALLOCATE(polyx(nx+4), polyz(nx+4), field(nfld), satnw2(nxyz), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
     PRINT *, 'ARRAY ALLOCATION FAILED: plotb2'
     ierr(187) = .TRUE.
     RETURN
  END IF
  !     ------------------------------------------------------------------
  IF(ABS(print_plotscalar) > 0._kdp) THEN
     CALL print_control(print_plotscalar,time,ltime,tchg,timprpscal,lprint)
     IF(lprint) THEN


        !         Output plotfile with the top rows first to match B2plot input
        !           write the header line for each time step for the pltout file
        IF (ioptpr(1) == 9) THEN
           WRITE(fupltsca) 'time=', year(time), nx, nz - 1, nfield, 1
        ELSE
           WRITE(fupltsca, '(A,G20.12,4I5)') 'time=', year(time), nx, nz - 1, nfield, 1
        END IF
        !       write the width of the region mesh, assuming left boundary is x=0
        width = xx(nx) + 0.5_kdp*dx(nx)
        xleft = -xx(1)
        xright = xx(nx) + dx(nx)
        ztop = zz(nz)
        ysl = IDNINT(14._kdp/11._kdp*akm(ztop))
        IF (ioptpr(1) == 9) THEN
           WRITE(fupltsca) akm(width), akm(xleft), akm(xright), ysl
        ELSE
           WRITE(fupltsca, '(4g20.12)') akm(width), akm(xleft), akm(xright), ysl
        END IF
        !        write labels for field variables and basin strata into interface
        IF (ioptpr(1) == 9) THEN
           WRITE(fupltsca) (ufield(i), ibdy(i), i=1, nfield)
        ELSE
           WRITE(fupltsca, '(4(a,i1))') (ufield(i), ibdy(i), i=1, nfield)
        END IF
        ubuff(1) = 'firstunit'
        IF (ioptpr(1) == 9) THEN
           WRITE(fupltsca) ubuff(1)
           WRITE(fupltsca) year(time)
        ELSE
           WRITE(fupltsca, '(a)') ubuff(1)
           WRITE(fupltsca, '(8g14.6)') year(time)
        END IF
        !             basin2 output for more than 1 horizon
        !       WRITE(fupltsca) (ubuff(n),n=1,nstrat)
        !       WRITE(fupltsca) (year(tdep(n)),n=1,nstrat)
        !         create a modified saturation array SATNW2 so that all supercritical
        !           nodes (P>Pc) have a saturation of -1 instead of 0.5
        DO  ic=1,nxyz
           IF (p(ic) > 220.55e6_kdp) THEN
              !              use the following if supercritical includes T>Tc
              !       IF ( (P(IC).GT.220.55D6.AND.SATNW(IC).EQ.0.5_KDP).OR.
              !    &       (TC(IC).GE.374.0_KDP.AND.SATNW(IC).EQ.0.0_KDP) ) THEN
              satnw2(ic) = -1._KDP
           ELSE
              satnw2(ic) = satnw(ic)
           END IF
        END DO
        !     Write data for each node to the plot dataset
        j = 1
        !                ----- top row, left boundary node -----
        i = 1
        ij = (j-1)*nx + i
        k = ktop(ij)
        ic = (k-1)*nx + i + (j-1)*nx*nz
        ik = (k-1)*nx + i
        ifield(1) = 0
        ifield(2) = 0
        ifield(3) = 2
        ifield(4) = (nz-k) + 1
        ifield(5) = -99
        ifield(6) = 1
        xposition = -akm(xx(1))
        zposition = ysl - akm(ztop) + akm(zz(k))
        field(1) = 0.0_kdp
        field(2) = 0.0_kdp
        field(3) = 0.0_kdp
        field(4) = akm(ztop-zz(k))
        field(5) = wpot(ic)
        field(6) = tc(ic)
        field(7) = atm(p(ic))
        field(8) = phi(ic)
        field(9) = LOG10(xk(ic)*1.01e8_kdp)
        field(10) = LOG10(zk(ic)*1.01e8_kdp)
        field(11) = satnw2(ic)
        field(12) = xkc(ic)*0.239e-7_kdp
        IF (ioptpr(1) == 9) THEN
           WRITE(fupltsca) (ifield(m), m=1, 6), xposition, zposition,  &
                (field(n), n=1, nfield)
        ELSE
           WRITE(fupltsca, 9005) (ifield(m), m=1, 6), xposition, zposition,  &
                (field(n), n=1, nfield)
9005       FORMAT (6I3, 8G14.6/(18X,8G14.6))
        END IF
        !        ----- top row, interior boundary nodes (top row of HYDROTHERM) -----
        DO  i = 1, nx
           ij = (j-1)*nx + i
           k = ktop(ij)
           ic = (k-1)*nx + i + (j-1)*nx*nz
           ik = (k-1)*nx + i
           ifield(1) = i
           ifield(2) = 0
           ifield(3) = 2
           ifield(4) = (nz-k) + 1
           ifield(5) = -99
           ifield(6) = 1
           xposition = akm(xx(i))
           zposition = ysl - akm(ztop) + akm(zz(k))
           field(1) = cmyr(xwvel(ic))
           field(2) = -cmyr(zwvel(ic))
           field(3) = 0.0_kdp
           field(4) = akm(ztop-zz(k))
           field(5) = wpot(ic)
           field(6) = tc(ic)
           field(7) = atm(p(ic))
           field(8) = phi(ic)
           field(9) = LOG10(xk(ic)*1.01e8_kdp)
           field(10) = LOG10(zk(ic)*1.01e8_kdp)
           field(11) = satnw2(ic)
           field(12) = xkc(ic)*0.239e-7_kdp
           IF (ioptpr(1) == 9) THEN
              WRITE(fupltsca) (ifield(m), m=1, 6), xposition, zposition,  &
                   (field(n), n=1, nfield)
           ELSE
              WRITE(fupltsca,9005) (ifield(m),m=1,6), xposition, zposition,  &
                   (field(n), n=1, nfield)
           END IF
        END DO
        !                ----- top row, right boundary node -----
        i = nx
        ij = (j-1)*nx + i
        k = ktop(ij)
        ic = (k-1)*nx + i + (j-1)*nx*nz
        ik = (k-1)*nx + i
        ifield(1) = nx + 1
        ifield(2) = 0
        ifield(3) = 2
        ifield(4) = (nz-k) + 1
        ifield(5) = -99
        ifield(6) = 1
        xposition = akm(xx(nx)+dx(nx))
        zposition = ysl - akm(ztop) + akm(zz(k))
        field(1) = 0.0_kdp
        field(2) = 0.0_kdp
        field(3) = 0.0_kdp
        field(4) = akm(ztop-zz(k))
        field(5) = wpot(ic)
        field(6) = tc(ic)
        field(7) = atm(p(ic))
        field(8) = phi(ic)
        field(9) = LOG10(xk(ic)*1.01e8_kdp)
        field(10) = LOG10(zk(ic)*1.01e8_kdp)
        field(11) = satnw2(ic)
        field(12) = xkc(ic)*0.239e-7_kdp
        IF (ioptpr(1) == 9) THEN
           WRITE(fupltsca) (ifield(m), m=1, 6), xposition, zposition,  &
                (field(n), n=1, nfield)
        ELSE
           WRITE(fupltsca, 9005) (ifield(m), m=1, 6), xposition, zposition,  &
                (field(n), n=1, nfield)
        END IF
        !              ########## main body of interior nodes ##########
        DO  k=nz-1,1,-1
           !                      --- left boundary node ---
           i = 1
           ij = (j-1)*nx + i
           ik = (k-1)*nx + i
           ikzp = ik + nx
           ifield(1) = 0
           ifield(2) = nz - k
           ifield(3) = 2
           IF (np(ik,j) == -1 .AND. k == ktop(ij)) ifield(3) = -1
           IF (np(ik,j) == 0) ifield(3) = -1
           ifield(4) = ifield(2) + 1
           IF (np(ik,j) <= 0 .AND. k >= ktop(ij)) ifield(4) = nz-ktop(ij)+1
           ifield(5) = ifield(2) - 1
           IF (np(ik,j) <= 0 .AND. k >= ktop(ij)) ifield(5) = 0
           IF (np(ikzp,j) <= 0 .AND. k+1 >= ktop(ij)) ifield(5) = 0
           ifield(6) = 1
           xposition = -1.0_KDP*akm(xx(1))
           zposition = ysl - akm(ztop) + akm(zz(k))
           IF (np(ik,j) <= 0 .AND. k >= ktop(ij))  &
                zposition = ysl - akm(ztop) + akm(zz(ktop(ij)-1))
           ic = (k-1)*nx + i + (j-1)*nx*nz
           IF (np(ik,j) == -1 .AND. k == ktop(ij)) ic = (k+1-1)*nx + i + (j-1)*nx*nz
           field(1) = 0.0_kdp
           field(2) = 0.0_kdp
           field(3) = 0.0_kdp
           field(4) = akm(ztop-zz(k))
           field(5) = wpot(ic)
           field(6) = tc(ic)
           field(7) = atm(p(ic))
           field(8) = phi(ic)
           field(9) = LOG10(xk(ic)*1.01e8_kdp)
           field(10) = LOG10(zk(ic)*1.01e8_kdp)
           field(11) = satnw2(ic)
           field(12) = xkc(ic)*0.239e-7_kdp
           IF (ioptpr(1) == 9) THEN
              WRITE(fupltsca) (ifield(m), m=1, 6), xposition, zposition,  &
                   (field(n), n=1, nfield)
           ELSE
              WRITE(fupltsca,9005) (ifield(m),m=1,6), xposition, zposition,  &
                   (field(n), n=1, nfield)
           END IF
           !        ----- interior nodes of each row (generated in HYDROTHERM) -----
           DO  i = 1, nx
              ij = (j-1)*nx + i
              ik = (k-1)*nx + i
              ikzp = ik + nx
              ifield(1) = i
              ifield(2) = nz - k
              ifield(3) = 0
              IF (np(ik,j) == -1 .AND. k /= ktop(ij)) ifield(3) = 2
              IF (np(ik,j) == -1 .AND. k == ktop(ij)) ifield(3) = -1
              IF (np(ik,j) == 0) ifield(3) = -1
              ifield(4) = ifield(2) + 1
              IF (np(ik,j) <= 0 .AND. k >= ktop(ij)) ifield(4)=nz-ktop(ij)+1
              ifield(5) = ifield(2) - 1
              IF (np(ik,j) <= 0 .AND. k >= ktop(ij)) ifield(5) = 0
              IF (np(ikzp,j) <= 0 .AND. k+1 >= ktop(ij)) ifield(5) = 0
              ifield(6) = 1
              xposition = akm(xx(i))
              zposition = ysl - akm(ztop) + akm(zz(k))
              IF (np(ik,j) <= 0 .AND. k >= ktop(ij))  &
                   zposition = ysl - akm(ztop) + akm(zz(ktop(ij)-1))
              ic = (k-1)*nx + i + (j-1)*nx*nz
              IF (np(ik,j) == -1 .AND. k == ktop(ij)) ic = (k+1-1)*nx + i + (j-1)*nx*nz
              field(1) = cmyr(xwvel(ic))
              field(2) = -cmyr(zwvel(ic))
              field(3) = 0.0_KDP
              field(4) = akm(ztop-zz(k))
              field(5) = wpot(ic)
              field(6) = tc(ic)
              field(7) = atm(p(ic))
              field(8) = phi(ic)
              field(9) = LOG10(xk(ic)*1.01e8_kdp)
              field(10) = LOG10(zk(ic)*1.01e8_kdp)
              field(11) = satnw2(ic)
              field(12) = xkc(ic)*0.239e-7_kdp
              IF (ioptpr(1) == 9) THEN
                 WRITE(fupltsca) (ifield(m), m=1, 6), xposition, zposition,  &
                      (field(n), n=1, nfield)
              ELSE
                 WRITE(fupltsca,9005) (ifield(m),m=1,6),xposition,zposition,  &
                      (field(n), n=1, nfield)
              END IF
           END DO
           !                        ----- right boundary node -----
           i = nx
           ij = (j-1)*nx + i
           ik = (k-1)*nx + i
           ikzp = ik + nx
           ifield(1) = nx + 1
           ifield(2) = nz - k
           ifield(3) = 2
           IF (np(ik,j) == -1 .AND. k == ktop(ij)) ifield(3) = -1
           IF (np(ik,j) == 0) ifield(3) = -1
           ifield(4) = ifield(2) + 1
           IF (np(ik,j) <= 0 .AND. k >= ktop(ij)) ifield(4) = nz-ktop(ij)+1
           ifield(5) = ifield(2) - 1
           IF (np(ik,j) <= 0 .AND. k >= ktop(ij)) ifield(5) = 0
           IF (np(ikzp,j) <= 0 .AND. k+1 >= ktop(ij)) ifield(5) = 0
           ifield(6) = 1
           xposition = akm(xx(nx)+dx(nx))
           zposition = ysl - akm(ztop) + akm(zz(k))
           IF (np(ik,j) <= 0 .AND. k >= ktop(ij))  &
                zposition = ysl - akm(ztop) + akm(zz(ktop(ij)-1))
           ic = (k-1)*nx + i + (j-1)*nx*nz
           IF (np(ik,j) == -1 .AND. k == ktop(ij)) ic = (k+1-1)*nx + i + (j-1)*nx*nz
           field(1) = 0.0_KDP
           field(2) = 0.0_KDP
           field(3) = 0.0_KDP
           field(4) = akm(ztop-zz(k))
           field(5) = wpot(ic)
           field(6) = tc(ic)
           field(7) = atm(p(ic))
           field(8) = phi(ic)
           field(9) = LOG10(xk(ic)*1.01E8_KDP)
           field(10) = LOG10(zk(ic)*1.01E8_KDP)
           field(11) = satnw2(ic)
           field(12) = xkc(ic)*0.239e-7_kdp
           IF (ioptpr(1) == 9) THEN
              WRITE(fupltsca) (ifield(m), m=1, 6), xposition, zposition,  &
                   (field(n), n=1, nfield)
           ELSE
              WRITE(fupltsca,9005) (ifield(m),m=1,6),xposition,zposition,  &
                   (field(n), n=1, nfield)
           END IF
        END DO
        !                  ########## bottom boundary row ##########
        k = 1
        !                      --- left boundary node ---
        i = 1
        ij = (j-1)*nx + i
        ik = (k-1)*nx + i
        ikzp = ik + nx
        ifield(1) = 0
        ifield(2) = nz
        ifield(3) = 2
        ifield(4) = -99
        ifield(5) = ifield(2) - 1
        IF (np(ik,j) <= 0 .AND. k >= ktop(ij)) ifield(5) = 0
        IF (np(ikzp,j) <= 0 .AND. k+1 >= ktop(ij)) ifield(5) = 0
        ifield(6) = 0
        xposition = -1.0_KDP*akm(xx(1))
        zposition = ysl - akm(ztop) + akm(zz(1)-dz(1))
        ic = (k-1)*nx + i + (j-1)*nx*nz
        field(1) = 0.0_KDP
        field(2) = 0.0_KDP
        field(3) = 0.0_KDP
        field(4) = akm(ztop-zz(1)+dz(1))
        field(5) = wpot(ic)
        !            linear extrapolation of temp below boundary
        field(6) = (tc(ic)-tc(ic+nx))* dz(1)/(0.5_KDP*dz(1)+0.5_KDP*dz(2)) + tc(ic)
        field(7) = atm(p(ic))
        field(8) = phi(ic)
        field(9) = LOG10(xk(ic)*1.01E8_KDP)
        field(10) = LOG10(zk(ic)*1.01E8_KDP)
        field(11) = satnw2(ic)
        field(12) = xkc(ic)*0.239e-7_kdp
        IF (ioptpr(1) == 9) THEN
           WRITE(fupltsca) (ifield(m), m=1, 6), xposition, zposition,  &
                (field(n), n=1, nfield)
        ELSE
           WRITE(fupltsca, 9005) (ifield(m), m=1, 6), xposition, zposition,  &
                (field(n), n=1, nfield)
        END IF
        !            ----- bottom boundary row, middle nodes -----
        DO  i = 1, nx
           ij = (j-1)*nx + i
           ik = (k-1)*nx + i
           ikzp = ik + nx
           ifield(1) = i
           ifield(2) = nz
           ifield(3) = 2
           ifield(4) = -99
           ifield(5) = ifield(2) - 1
           IF (np(ik,j) <= 0 .AND. k >= ktop(ij)) ifield(5) = 0
           IF (np(ikzp,j) <= 0 .AND. k+1 >= ktop(ij)) ifield(5) = 0
           ifield(6) = 0
           xposition = akm(xx(i))
           zposition = ysl - akm(ztop) + akm(zz(1)-dz(1))
           ic = (k-1)*nx + i + (j-1)*nx*nz
           field(1) = 0.0_KDP
           field(2) = 0.0_KDP
           field(3) = 0.0_KDP
           field(4) = akm(ztop-zz(1)+dz(1))
           field(5) = wpot(ic)
           !            linear extrapolation of temp below boundary
           field(6) = (tc(ic)-tc(ic+nx))* dz(1)/(0.5_KDP*dz(1)+0.5_KDP*dz(2)) + tc(ic)
           field(7) = atm(p(ic))
           field(8) = phi(ic)
           field(9) = LOG10(xk(ic)*1.01E8_KDP)
           field(10) = LOG10(zk(ic)*1.01E8_KDP)
           field(11) = satnw2(ic)
           field(12) = xkc(ic)*0.239e-7_kdp
           IF (ioptpr(1) == 9) THEN
              WRITE(fupltsca) (ifield(m), m=1, 6), xposition, zposition,  &
                   (field(n), n=1, nfield)
           ELSE
              WRITE(fupltsca,9005) (ifield(m),m=1,6),xposition,zposition,  &
                   (field(n), n=1, nfield)
           END IF
        END DO
        !            ----- bottom boundary row, right node -----
        i = nx
        ij = (j-1)*nx + i
        ik = (k-1)*nx + i
        ikzp = ik + nx
        ifield(1) = nx + 1
        ifield(2) = nz
        ifield(3) = 2
        ifield(4) = -99
        ifield(5) = ifield(2) - 1
        IF (np(ik,j) <= 0 .AND. k >= ktop(ij)) ifield(5) = 0
        IF (np(ikzp,j) <= 0 .AND. k+1 >= ktop(ij)) ifield(5) = 0
        ifield(6) = 0
        xposition = akm(xx(nx)+dx(nx))
        zposition = ysl - akm(ztop) + akm(zz(1)-dz(1))
        ic = (k-1)*nx + i + (j-1)*nx*nz
        field(1) = 0.0_KDP
        field(2) = 0.0_KDP
        field(3) = 0.0_KDP
        field(4) = akm(ztop-zz(1)+dz(1))
        field(5) = wpot(ic)
        !            linear extrapolation of temp below boundary
        field(6) = (tc(ic)-tc(ic+nx))* dz(1)/(0.5_KDP*dz(1)+0.5_KDP*dz(2)) + tc(ic)
        field(7) = atm(p(ic))
        field(8) = phi(ic)
        field(9) = LOG10(xk(ic)*1.01E8_KDP)
        field(10) = LOG10(zk(ic)*1.01E8_KDP)
        field(11) = satnw2(ic)
        field(12) = xkc(ic)*0.239e-7_kdp
        IF (ioptpr(1) == 9) THEN
           WRITE(fupltsca) (ifield(m), m=1, 6), xposition, zposition,  &
                (field(n), n=1, nfield)
        ELSE
           WRITE(fupltsca, 9005) (ifield(m), m=1, 6), xposition, zposition,  &
                (field(n), n=1, nfield)
        END IF
        !       ----------------------------------------------------------
        !     write bounding polygon.  the indexing for polygon coordinates
        !      (xline,zline) begins at the upper left, then lower left,
        !      continues along the base of the cross-section, lower right,
        !      upper right, and returns along the basin surface.
        npoly = 2*nx + 4
        polyx(1) = 0.0_kdp
        polyz(1) = akm(zztop(1)-ztop) + ysl
        polyx(2) = 0.0_kdp
        polyz(2) = ysl - akm(ztop)
        DO  i = 1, nx
           polyx(i+2) = akm(xx(i))
           polyz(i+2) = ysl - akm(ztop)
           polyx(npoly+1-i) = akm(xx(i))
           polyz(npoly+1-i) = akm(zztop(i)-ztop) + ysl
        END DO
        polyx(nx+3) = akm(xx(nx)+0.5_KDP*dx(nx))
        polyz(nx+3) = ysl - akm(ztop)
        polyx(nx+4) = akm(xx(nx)+0.5_KDP*dx(nx))
        polyz(nx+4) = akm(zztop(nx)-ztop) + ysl

        IF (ioptpr(1) == 9) THEN
           WRITE(fupltsca) (polyx(n), polyz(n), n=1, npoly)
        ELSE
           WRITE(fupltsca, '(6f15.7)') (polyx(n), polyz(n), n=1, npoly)
        END IF
        !         write positions of contacts between formations.  scan from
        !                oldest (basal) strata upward.
        IF (ioptpr(1) == 9) THEN
           WRITE(fupltsca) 0
        ELSE
           WRITE(fupltsca, '(i5)') 0
        END IF
     END IF
  END IF
  DEALLOCATE(polyx, polyz, field, satnw2)
END SUBROUTINE plotb2
