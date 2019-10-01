SUBROUTINE init_ht(flag)
  ! ... Initializes the simulation run and executes the program up to the
  ! ...      time stepping loop
  USE bc
  USE control
  USE fdeq
  USE mesh
  USE parameters, ONLY: lrxftreq, lrxftimreq, pwb
  USE variables
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: flag
  !
  CHARACTER(LEN=80) :: ufile
  INTEGER :: ic, ie
  INTEGER :: i, ik, ikl, iklxp, ikzp, j, k, mxm, mym, mzm
  INTEGER :: m
  REAL(KIND=kdp) :: xmsum, ensum
  INTERFACE
     SUBROUTINE phreg(p,h)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
     END SUBROUTINE phreg
     SUBROUTINE tempdegc(p,h,tc)
       USE machine_constants, ONLY: kdp
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: p, h
       REAL(KIND=kdp), DIMENSION(:), INTENT(INOUT) :: tc
     END SUBROUTINE tempdegc
  END INTERFACE
  ! ...  Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.2 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  ! ... initialize variables
!!$  mrno = 0
  resmold = 0._kdp
  resmold2 = 0._kdp
  reseold = 0._kdp
  reseold2 = 0._kdp
  resmaxmall = 0._kdp
  ltimresmaxm = 0
  resmaxeall = 0._kdp
  iconvrg = 0
  iphschgtot = 0
  ierr = .FALSE.
  flag = 0
  kod = 0
  ltimresmaxe = 0
  lnri = 0
  lnritotl = 0
  !$$  slmeth = 0
  tot_iter_gm = 0
  lrxftimreq = .FALSE.
  lrxftreq = .FALSE.
  ilnrimax = 0
  iter_gm_max = 0
  ltime = 0
  ltmcntl = 0
  itmcntl = 0
  unconfined = .FALSE.
  pwb = 0._kdp
  wrelaxp = 1._kdp
  wrelaxh = 1._kdp
  time = 0._kdp
  tchg = 1.e99_kdp     ! ... To set print stops initially
  ioptpr = 0
  totfi = 0._kdp
  totfp = 0._kdp
  tothi = 0._kdp
  tothp = 0._kdp
  tfres = 0._kdp
  thres = 0._kdp
  tcfsbc = 0._kdp
  tcfpfbc = 0._kdp
  tcfsepbc = 0._kdp
  totwfi = 0._kdp
  totwfp = 0._kdp
  tchsbc = 0._kdp
  tchpfbc = 0._kdp
  tchsepbc = 0._kdp
  tchhcbc = 0._kdp
  totwhi = 0._kdp
  totwhp =0._kdp
  lpltexpl = 0
  ! ... Read static data
  CALL gdata(flag)
  DO ie=1,200
     IF(ierr(ie)) THEN
        flag = 1
     END IF
  END DO
  IF(flag > 0) RETURN
  ! ... Save the maximum values of initial pressure and initial enthalpy fields
  maxpressic = MAXVAL(p)
  maxenthic = MAXVAL(h)
  CALL tempdegc(p,h,tc)     ! ... Calculate temperature for i.c. P and H
  nppn = 0          ! ... counter for number of times SINK has been called
  CALL sink(nppn,flag)     ! ... read initial source/sink data
  IF(flag > 0) RETURN     ! ... exit if a read or i.c. problem in sink
  maxdelt = tchg     ! ... *** temporary 
  ltmcntl = 6
  CALL tcalc  ! ... Calculate transmissivity and thermal conductivity for cell faces
  CALL phreg(p,h)  ! ... Determine initial phase region (water, 2 phase, steam, or air-water) for cells
  CALL enthrock  ! ... Calculate enthalpy of the porous matrix
  CALL prpty(0)  ! ... Calculate fluid property coefficients and derivatives
  CALL wellallo  ! ... Allocate the heat and mass of source/sink terms to the cells
  ! ... Initialize spatial weight factors for midpoint weighting
  DO  mxm = 1, nxx*ny*nz
     usx(mxm) = ax(mxm)
     uwx(mxm) = ax(mxm)
  END DO
  DO  mzm = 1, nx*ny*nzz
     usz(mzm) = az(mzm)
     uwz(mzm) = az(mzm)
  END DO
  IF (ny > 1) THEN
     DO  mym = 1, nx*nyy*nz
        usy(mym) = ay(mym)
        uwy(mym) = ay(mym)
     END DO
  END IF
  IF (ioptupst) CALL upstre     ! ... Determine the upstream node for each cell
  CALL storativ          ! ... Calculate the initial mass and energy storage terms (capacitances)
  ! ... Save mass and energy amounts from the initial conditions (arrays)
  xmoldt = xm
  enoldt = en
  ! ...    sum the total mass and energy in the problem domain
  xmsum = 0._kdp
  ensum = 0._kdp
  DO  ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF (ibc(ik,j) /= -1) THEN        ! ... Select active cells
        xmsum = xmsum + xm(ic)
        ensum = ensum + en(ic)
     END IF
  END DO
  fir0 = xmsum
  ehir0 = ensum
  ! ... Write out the initial conditions
  ! ... Set the print masks for interfacial array printouts
  !                  ------------ X-direction ------------
  lprntxif = -1      ! ... Initialize all interfaces out of the active region
  DO  ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ! ... indices in x-z slice
     ik = (k-1)*nx + i
     ikl = (k-1)*(nx+1) + i
     iklxp = ikl + 1
     IF (ibc(ik,j) /= -1) THEN
        ! ... Set print mask to print facial values for both sides of each cell
        ! ...      included in the region
        lprntxif(ikl,j) = 0
        lprntxif(iklxp,j) = 0
     END IF
  END DO
  !                  ------------ Y-direction ------------
  lprntyif = -1      ! ... Initialize all interfaces out of the active region
  DO  ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ! ... indices in x-z slice
     ik = (k-1)*nx + i
     IF (ibc(ik,j) /= -1) THEN
        ! ... Set print mask to print facial values for both sides of each cell
        ! ...      included in the region
        lprntyif(ik,j) = 0
        lprntyif(ik,j+1) = 0
     END IF
  END DO
  !                  ------------ Z-direction ------------
  lprntzif = -1      ! ... Initialize all interfaces out of the active region
  DO  ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ! ... indices in x-z slice
     ik = (k-1)*nx + i
     ikzp = ik + nx
     IF (ibc(ik,j) /= -1) THEN
        ! ... Set print mask to print facial values for both sides of each cell
        ! ...      included in the region
        lprntzif(ik,j) = 0
        lprntzif(ikzp,j) = 0
     END IF
  END DO
  IF (ABS(print_vel) > 0._kdp .OR. ABS(print_plotvector) > 0._kdp .OR.  &
       ABS(print_plotscalar) > 0._kdp .OR.  &
       ioptprta(9) >= 1 .OR. ioptprta(10) >= 1 .OR.  &
       ioptprta(11) >= 1 .OR. ioptprta(15) >= 1 .OR.  &
       ioptprta(16) >= 1 .OR. ioptprta(17) >= 1 .OR. ioptprta(27) >=1) THEN
     ! ... Calculate the water and steam velocity for printing if necessary
     CALL velocity
  END IF
  CALL pdata(flag)     ! ... Write tabular information to output files
  IF(flag > 0) RETURN
  SELECT CASE (ioptpr(1))
  CASE (1)        ! ... write Explorer plotfile
     ufile = '100'
     CALL gfiles(6, ufile, flag)
     IF(flag > 0) RETURN
     CALL plotexpl
  CASE (2,3)      ! ... write gnuplot plotfile
     ufile = '100'
     CALL gfiles(7, ufile, flag)
     IF(flag > 0) RETURN
     CALL plotgnu
  CASE (4,5)      ! ... write HTView IDL plotfile
     ufile = '100'
     CALL gfiles(7, ufile, flag)
     IF(flag > 0) RETURN
     CALL plotidl
  CASE (6)        ! ... Write columnar xyz plotfile
     ufile = 'xyz'
     CALL gfiles(7, ufile, flag)
     IF(flag > 0) RETURN
     CALL plotxyz
  CASE (8,9)      ! ... write B2plotfile
     ufile = '100'
     CALL gfiles(7, ufile, flag)
     IF(flag > 0) RETURN
     CALL plotb2
  END SELECT
  ! ... Set the time for next printout
  prtimenxt = MIN(tchg, timprp, timprh, timprt, timprsat, timprd, timprv, timprpot,  &
       timprvel, timprbcf, timprsrc, timprprop, timprpor, timprperm, timprbal, timprdump,  &
       timprresid, timprdimno, timprpscal, timprpvec, timprtemp)
  iplot = 0
  ! ... If close to time for printout, move to time for printout
  ! ... If overshot time for printout, back up
  IF(ABS(time+delt-prtimenxt) <= 0.1*delt .OR. time+delt > prtimenxt) THEN
     delt = prtimenxt - time
     iplot = 1
     ltmcntl = 8
  END IF
!  IF (prfreq < 0._kdp) iplot = 1
  icall = 0
  ! ... If close to time for change, move to time for change or
  ! ...     if overshot time for change, back up
  IF(ABS(time+delt-tchg) <= 0.2*delt .OR. time+delt > tchg) THEN
     delt = tchg - time
     icall = 1
     ltmcntl = 7
  END IF
  ! ... Save the phase indicator and time step?? from the initial conditions
  DO  ic = 1,nxyz
     indoldnr(ic) = ind(ic)
  END DO
  deltoldt = delt
  deltmin = delt
  deltmax = delt
  ! ... Check for error flags
  DO  ie=1,200
     IF(ierr(ie)) flag = 1
  END DO
END SUBROUTINE init_ht
