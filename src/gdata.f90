SUBROUTINE gdata(flag)
  ! ... Purpose:  To read and write problem definition information, set up the
  ! ...      initial conditions, and read and write simulation control information
  ! ...      Allocates the array space needed for this simulation
  USE machine_constants, ONLY: kdp
  USE f_units
  USE bc
  USE parameters          !!$*** for flint
  USE control
  USE fdeq
  USE i_c
  USE ilupc_mod
  USE mesh
!!$  USE parameters
  USE solver
  USE solver_gmres
  USE variables
  USE units
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: flag
  !
!!$  CHARACTER(LEN=10), EXTERNAL :: uppercase
  CHARACTER(LEN=80) :: upltdir, ufile, ustring
  CHARACTER(LEN=20) :: date_time
  !
  CHARACTER(LEN=80) :: cfmta, cfmtb, cfmtc, cfmtd
  CHARACTER(LEN=22), PARAMETER :: cfmta1="(/(8(1pg9.3,tr1)))   ",  &
       cfmta2="(/(10(1pg9.3,tr1)))  ",  &
       cfmta3="(/(16(1pg9.3,tr1)))  ",  &
       cfmtb1="(/tr7,A/(8(i5,tr5))) ",  &
       cfmtb2="(/tr7,A/(10(i5,tr5)))",  &
       cfmtb3="(/tr7,A/(16(i5,tr5)))"
  CHARACTER(LEN=80), PARAMETER :: cfmtc1="(/,'Column (I)->',i3,13I5/(12X,13I5))  ",  &
       cfmtc2="(/'Column (I)->',i3,23I5/(12X,23I5))    ",  &
       cfmtc3="(/'Column (I)->',i3,28I5/(12X,28I5))    ",  &
       cfmtd1="('Row',i3,' -> ',14I5/(12X,13I5))       ",  &
       cfmtd2="('Row',i3,' -> ',24I5/(12X,23I5))       ",  &
       cfmtd3="('Row',i3,' -> ',29I5/(12X,28I5))       "
  INTEGER :: i, ibctx, ic, icdum, icxf, icxl, icxzf, icxzl, ie, ii,  &
       ij, ik, iknext, ikxf, ikxl, inext, irx
  INTEGER :: j, k, knext, kodparm, ldum,  &
       lll, lrowcol, ml, mlnext, narray, nby,  &
       nconst, ndumy, nodseqno, nodseqnoglb, nout
  ! ...   nodseqno - Node sequence number counter for one x-z slice
  ! ...   nodseqnoglb - Global node sequence number counter for 3d iterative solver
  INTEGER :: a_err
  INTEGER :: i1, mprtp1
  REAL(KIND=kdp) :: theta
  REAL(KIND=kdp) :: gveqnm, gveqne, gveqnead, gveqnecd, qfluxsp, qhfluxsp
!!$  INTEGER, EXTERNAL :: iargc
  LOGICAL :: newfmt, lyr
  INTERFACE
     FUNCTION uppercase(string) RESULT(outstring)
       CHARACTER(LEN=*), INTENT(IN) :: string
       CHARACTER(LEN=LEN(string)) :: outstring
     END FUNCTION uppercase
  END INTERFACE     !!$*** for flint
  INTERFACE         !!$*** for flint
     SUBROUTINE printar(ndim,array,lprnt,fu,cnv,jfmt,nnoppr)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: ndim
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: array
       INTEGER, DIMENSION(:,:), INTENT(IN) :: lprnt
       INTEGER, INTENT(IN) :: fu
       REAL(KIND=kdp), INTENT(IN) :: cnv
       INTEGER, INTENT(IN) :: jfmt
       INTEGER, INTENT(IN) :: nnoppr
     END SUBROUTINE printar
  END INTERFACE
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.5 $//$Date: 2010/03/30 23:55:48 $'
  ! ... ------------------------------------------------------------------
  ! ... open main input and output files
  ufile = ' '
  CALL gfiles(1,ufile,flag)
  ! ... Strip out comment lines from input file
  CALL rcomment
  WRITE(fupdef,'(A)')  &
       '***** HYDROTHERM - Water, Steam, Two-phase, and Supercritical Flow Simulator *****'
  READ(fuins,'(A)') title1
  READ(fuins,'(A)') title2
  ! ... check first word of input file for "TITLE"
  ustring = ADJUSTL(title1)
  ustring = uppercase(ustring(1:5))
  IF (ustring(1:5) /= 'TITLE') THEN
     WRITE(fustdout,'(/a)') '*** ERROR: Invalid First Line of Input File ***'
     WRITE(fuclog,'(/a)') '*** ERROR: Invalid First Line of Input File ***'
     ierr(6) = .TRUE.
  END IF
  ! ... Identify the format of input file
  ustring = ADJUSTL(title2)
  ! ... Extract the version name
  i1 = INDEX(ustring,' ')
  ustring = uppercase(ustring(1:i1-1))
  informtp = 0
  SELECT CASE (ustring(1:i1-1))
  CASE ("2FOR")
     informt = 1
  CASE ("VER2")
     informt = 2
  CASE ("VER3")
     informt = 3
  CASE ("VER3.2")
     informt = 3
     informtp = 2
  CASE default
     WRITE(fustdout,9005) 'ERROR: *** INVALID FORMAT VERSION OF INPUT FILE ***',  &
          'Line 2 begins with: ',title2(1:10)
     WRITE(fuclog,9005) 'ERROR: *** INVALID FORMAT VERSION OF INPUT FILE ***',  &
          'Line 2 begins with: ',title2(1:10)
9005 FORMAT(tr5,a/tr10,a,a10)
     ierr(7) = .TRUE.
     RETURN
  END SELECT
  ! ...  Extract title string 
  title1 = TRIM(title1(6:))
  title2 = TRIM(title2(i1:))
  WRITE(fupdef,'(2X,76(''-'')/2X,A/2X,A/2X,76(''-''))') title1, title2
  ! ...         program and run information 
  WRITE(fupdef,9235) '--- Run Information ---'
  WRITE(fuclog,9235) '--- Run Information ---'
  CALL timestamp(date_time)
  WRITE(fupdef,9010) 'Execution Date and Time:   ',date_time
9010 FORMAT (tr10,2a)
  WRITE(fupdef,9180) 'Input filename:  ',TRIM(ufile)
  WRITE(fuclog,9180) 'Input filename:  ',TRIM(ufile)
  WRITE(fupdef,9175) 'Input format version:  '//ustring(1:8)
  WRITE(fuclog,9175) 'Input format version:  '//ustring(1:8)
  WRITE(fupdef,9180) 'Output filename suffix for this run:  ',usuff
  WRITE(fuclog,9180) 'Output filename suffix for this run:  ',usuff
  WRITE(fuclog,9010) 'Execution Date and Time:   ',date_time
  ! ... write the name of the directory for the explorer plotfiles
  ! ... set up by the third argument from the command line
  ! ... Temporarily inactivated
  upltdir = ' '
!!$  !..      narg = IARGC()
!!$  !..      IF (narg.EQ.3) THEN
!!$  !..        upltdir = ' '
!!$  !..        CALL GETARG(3, upltdir)
!!$  !..        CALL GETCWD(dirname)
!!$  !..        defltdir =
!!$  !..     &  TRIM(dirname)//'/'//TRIM(upltdir)//'/'
!!$  !..        WRITE(fuclog,9180) 'Explorer plotfiles directory: ',
!!$  !..     &                  TRIM(defltdir)
!!$  !..      ENDIF
  WRITE(fupdef,9235) '--- Program Version ---'
  WRITE(fupdef,9335) 'HYDROTHERM  Version:'//version_name
  !$$       'Release date:  unreleased'
  WRITE(fuclog,9235) '--- Program Version ---'
  WRITE(fuclog,9335) 'HYDROTHERM  Version:'//version_name
  !$$       'Release date:  unreleased'
  WRITE(fustdout,9335) 'HYDROTHERM  Version:'//version_name
9335 FORMAT(TR10,A/TR10,A)
  ! ...         ----- dimensions of problem domain -----
  ! ...   Nx - number of cell columns (X-dimension)
  ! ...   Ny - number of cell vertical slices (Y-dimension)
  ! ...   Nz - number of cell rows (Z-dimension)
  ! ...   Mbw - bandwidth, if set to 0, program will calculate
  ! ...             maximum value (ignores specified value nodes)
  ! ...   Timestr - starting time for this run.  Non-zero used to restart old runs
  ! ...   Lyr(logical)  - F for seconds; T for years
  ! ...   Iyr (integer) - 1 for time input in seconds; 2 for time input in years
  mbw = 0
  SELECT CASE (informt)
  CASE (1)
     READ(fuins,9000)
     READ(fuins,*) nx, ny, nz, mbw, timestr, lyr
     iyr=1
     IF(lyr) iyr=2
  CASE (2)
     READ(fuins,*,ERR=390) nx, ny, nz, timestr, lyr
     mbw = 0
     iyr=1
     IF(lyr) iyr=2
  CASE (3)
     READ(fuins,*,ERR=390) nx, ny, nz, timestr, iyr
  END SELECT
  ! ... Convert start time in years to start time in seconds (cgs units)
  IF(iyr == 2) timestr = sec(timestr)
  time = timestr
  WRITE(fupdef,9015) '--- Dimensions of Problem Domain and Starting Time ---',  &
       'Input Value','Number of:  Cell Columns (X) ..... ',nx,  &
       '            Cell Rows (Z) ........ ',nz,  &
       '            Cell Slices (Y) ...... ',ny
9015 FORMAT (/A/TR31,A/TR10,A,I5/TR10,A,I5/TR10,A,I5)
  CALL alloc_mesh(a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: gdata, alloc_mesh"
     ierr(199) = .TRUE.
     RETURN
  ENDIF
  ! ... check for dimension errors
  CALL alloc_fdeq(a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: gdata, alloc_fdeq"
     ierr(199) = .TRUE.
     RETURN
  ENDIF
  IF (iyr == 2) THEN
     WRITE(fupdef,9020) 'Initial Time .......... ',year(time),' (yr)'
     WRITE(fuclog,9020) 'Initial Time .......... ',year(time),' (yr)'
9020 FORMAT (tr10,a,1pg10.4,a)
     WRITE(fupdef,9240) 'All Time Variables Input in Years'
9240 FORMAT (tr10,a)
  ELSE
     WRITE(fupdef,9020) 'Initial Time .......... ',time,' (s)'
     WRITE(fuclog,9020) 'Initial Time .......... ',time,' (s)'
     WRITE(fupdef,9240) 'All Time Variables Input in Seconds'
  END IF
  nxx = nx + 1
  nyy = ny + 1
  nzz = nz + 1
  ! ...           ----- time stepping control -----
  ! ...   Ltimemax - maximum number of time steps
  ! ...   Tmincmx - maximum factor by which time step may be increased
  ! ...   Pchgmxts - max allowable % change in pressure for a node per time step; confined flow
  ! ...            - max allowable change in pressure for a node per time step; unconfined flow
  ! ...   Hchgmxts - max allowable % change in enthalpy for a node per time step
  ! ...   Wschgmx - max absolute change in water sat'n for a node per time step
  ! ...   Flwrlx - max number of cells through which fluid may flow per time step
  ! ...   mindelt - minimum time step length allowed
  ! ...   maxdelt - maximum time step length allowed
  ! ...   Ltmcutmx - number of times the delta time may be cut during a given time step
  IF (informt == 1) READ(fuins,9000)
  READ(fuins,*) ltimemax, tmincmx, pchgmxts, hchgmxts, wschgmx, flwrlx, mindelt, ltmcutmx
!!$  READ(fuins,*) ltimemax, tmincmx, pchgmxts, hchgmxts, wschgmx, flwrlx,  &
!!$        mindelt, maxdelt, ltmcutmx  *****activate later
  IF(iyr == 2) mindelt = sec(mindelt)
!!$  IF(iyr == 2) maxdelt = sec(maxdelt)
  ! ... Set defaults if necessary
  IF(tmincmx <= 1._kdp) tmincmx = 1.3_kdp
  !*****unconfined is not defined yet so this logic can not work. Later will be replaced
  ! ... ****** by other criteria independent of confinement. *****
  IF(unconfined .AND. pchgmxts <= 0._kdp) pchgmxts = 1.e5_kdp
  IF(.NOT.unconfined .AND. pchgmxts <= 0._kdp) pchgmxts = 10.0_kdp  ! ... % change
  IF(hchgmxts <= 0._kdp) hchgmxts = 5.0_kdp   ! ... % change
  IF(wschgmx <= 0._kdp) wschgmx = 0.03_kdp
  flwrlx = 1.0_kdp     ! ... forced value for now
  IF(ltmcutmx == 0) ltmcutmx = 10
  ! ***** unconfined is always false here, not read yet***
  IF(unconfined) THEN
     WRITE(fupdef,9125) ltimemax, tmincmx, pchgmxts, hchgmxts, wschgmx
9125 FORMAT (/, '--- Time Step Control ---',/TR10,  &
          'Maximum number of time steps:                    ',i7/tr10,  &
          'Max factor for increasing time step:               ',f5.2/tr10,  &
          'Max change in nodal pressure per time step:        ',1pg10.3/tr10,  &
          'Max % change in nodal enthalpy per time step:      ',f5.1/tr10,  &
          'Max change in nodal water satn per time step:      ',f5.3)
  ELSE
     WRITE(fupdef,9025) ltimemax, tmincmx, pchgmxts, hchgmxts, wschgmx
9025 FORMAT (/, '--- Time Step Control ---', /,TR10,  &
          'Maximum number of time steps:                    ',i7/TR10,  &
          'Max factor for increasing time step:               ',f5.2/TR10,  &
          'Max % change in nodal pressure per time step:      ',f5.1/TR10,  &
          'Max % change in nodal enthalpy per time step:      ',f5.1/TR10,  &
          'Max change in nodal water satn per time step:      ',f5.3)
  END IF
  IF (iyr == 2) THEN
     WRITE(fupdef,9120) 'Minimum Time Step Limit:        ',year(mindelt),' (yr)'
!!$  WRITE(fupdef,9120) 'Maximum Time Step Limit:        ',year(maxdelt),' (yr)'
  ELSE
     WRITE(fupdef,9120) 'Minimum Time Step Limit:        ',mindelt,' (s)'
!!$  WRITE(fupdef,9120) 'Maximum Time Step Limit:        ',maxdelt,' (s)'
9120 FORMAT (tr10,a,tr18,1pg10.4,a)
  END IF
  WRITE(fupdef,9030)  &
       'Max number of time step cuts per time step:        ',ltmcutmx
9030 FORMAT (tr10,a,i5)
  tmincmx = MAX(tmincmx,1._kdp)
  ! ...                ----- Newton-Raphson iteration controls -----
  ! ...    Lnrimax - maximum number of newton-raphson iterations
  ! ...    Resmas - max residual mass flux (g/s-cm^3) allowable for N-R
  ! ...      convergence
  ! ...    Reseng - max residual energy flux (erg/s-cm^3) allowable for N-R
  ! ...      convergence
  ! ... Iphschgmx - max number of nodes allowed to change phase per iteration
  ! ... Pchgmxnr - max % P change between NR iter for convergence 2ndary crit.; confined flow
  ! ...          - max unscaled P change between NR iter for convergence 
  ! ...               2ndary crit.; unconfined flow
  ! ... Hchgmxnr - max % H change between NR iter for convergence 2ndary crit.
  IF (informt == 1) THEN
     READ(fuins,9000)
     READ(fuins,*) lnrimax, resmas, reseng, iphschgmx
     pchgmxnr = 0.1_kdp     ! ... Forced values for old input format
     hchgmxnr = 0.1_kdp
  ELSE
     READ(fuins,*) lnrimax, resmas, reseng, iphschgmx, pchgmxnr, hchgmxnr
  END IF
  ! ... Set defaults if necessary
  IF(lnrimax == 0) lnrimax = 10
  IF(resmas <= 0._kdp) resmas = 1.e-12_kdp 
  IF(reseng <= 0._kdp) reseng = 1.e-2_kdp
  IF(iphschgmx == 0) iphschgmx = 10
  IF(unconfined .AND. pchgmxnr <= 0._kdp) pchgmxnr = 1.e5_kdp   ! ... 0.001*P at 1 km
  IF(.NOT.unconfined .AND. pchgmxnr <= 0._kdp) pchgmxnr = 0.1_kdp   ! ... % change
  IF(hchgmxnr <= 0._kdp) hchgmxnr = 0.1_kdp   ! .... % change
  IF(unconfined) THEN
     WRITE(fupdef,9135) lnrimax, resmas, reseng, iphschgmx, pchgmxnr, hchgmxnr
9135 FORMAT (/, '--- Newton-Raphson Iteration Convergence Criteria ---'  &
          ,/, 1P,TR10, 'Max number of iterations per time step:            ',i5/  &
          TR10, 'Max residual mass change (g/s-cm^3):             ',e11.2/  &
          TR10, 'Max residual energy change (erg/s-cm^3):         ',e11.2/  &
          TR10, 'Max number of nodes allowed to change phase per NR '/  &
          tr15,'iteration:                                    ',i5/ 5X,  &
          '-- Secondary Criteria (apply when residual criteria fail) --'/TR10,  &
          'Max unscaled pressure change per NR iteration:          ',1pe10.3/TR10,  &
          'Max % enthalpy change per NR iteration:          ',0PF5.1)
  ELSE
     WRITE(fupdef,9035) lnrimax, resmas, reseng, iphschgmx, pchgmxnr, hchgmxnr
9035 FORMAT (/, '--- Newton-Raphson Iteration Convergence Criteria ---'  &
          ,/, 1P,TR10, 'Max number of iterations per time step:            ',i5/  &
          TR10, 'Max residual mass change (g/s-cm^3):             ',e11.2/  &
          TR10, 'Max residual energy change (erg/s-cm^3):         ',e11.2/  &
          TR10, 'Max number of nodes allowed to change phase per NR '/  &
          tr15,'iteration:                                    ',i5/ 5X,  &
          '-- Secondary Criteria (apply when residual controls fail) --'/TR10,  &
          'Max % pressure change per NR iteration:          ',0PF5.1/TR10,  &
          'Max % enthalpy change per NR iteration:          ',0PF5.1)
  END IF
  ! ... Select the linear equation solver
  IF (informt == 3) THEN
     READ(fuins,*) slmeth
  ELSE
     slmeth = 1
  END IF
  ! ... check for errors
  IF (slmeth < 1 .OR. slmeth > 2) THEN          ! ... Unknown solver selection
     WRITE(fustdout,9250) '*** ERROR - Invalid index for linear equation solver  ***'
     WRITE(fuclog,9250) '*** ERROR - Invalid index for linear equation solver  ***'
     ierr(128) = .TRUE.
  END IF
  IF(((nx == 1 .AND. ny == 1) .OR. (ny == 1 .AND. nz == 1)) .AND. slmeth /= 1) THEN
     WRITE(fustdout,9250) '*** WARNING - One-dimensional region, '//  &
          'direct solver will be used  ***'
     WRITE(fuclog,9250) '*** WARNING -  One-dimensional region, '//  &
          'direct solver will be used  ***'
     slmeth = 1
!!$ ... ierr(129) = .TRUE.
  END IF
!  IF((nx > 1 .AND. ny > 1 .AND. nz > 1) .AND. slmeth == 1) THEN
!     WRITE(fustdout,9250) '*** ERROR - SSOR solver no longer supported for '//  &
!          '3-dimensional regions; Use GMRES solver  ***'
!     WRITE(fuclog,9250)  '*** ERROR - SSOR solver no longer supported for '//  &
!          '3-dimensional regions; Use GMRES solver  ***'
!     ierr(129) = .TRUE.
!  END IF
  IF((nx > 1 .AND. ny > 1 .AND. nz == 1) .AND. slmeth == 1) THEN
     WRITE(fustdout,9250) '*** ERROR - SSOR solver no longer supported for '//  &
          '2-dimensional areal regions; Use GMRES solver  ***'
     WRITE(fuclog,9250)  '*** ERROR - SSOR solver no longer supported for '//  &
          '2-dimensional areal regions; Use GMRES solver  ***'
     ierr(129) = .TRUE.
  END IF
  WRITE(fuclog,'(a)') '--- Linear Equation Solver ---'
  SELECT CASE (slmeth)
  CASE (1)
     WRITE(fupdef,9235) '--- SSOR Linear Equation Solver selected'        
     WRITE(fuclog,9235) '--- SSOR Linear Equation Solver selected'        
     ! ...        ----- Slice-Successive Over Relaxation Control -----
     ! ...    Lssormax - maximum number of SSOR iterations / timestep
     ! ...    Tol - tolerance (usually 0.01)
     ! ...    Wo - ssor relaxation factor (usually 1)
     IF (informt == 1) READ(fuins,9000)
     READ(fuins,*) lssormax, tol, wo
     ! ... Set defaults if necessary
     IF(lssormax <= 0) lssormax = 100
     IF(tol <= 0._kdp) tol = 1.e-5_kdp
     IF(wo <= 0._kdp) wo = 1.3_kdp
     WRITE(fupdef,9235) '--- Slice-Successive Over Relaxation Solver Control ---'
     IF (ny == 1) THEN
        WRITE(fupdef,9175) 'Direct linear solver used for 2-d slices; '//  &
             'No SSOR iterations necessary'
     ELSE
        WRITE(fupdef,9040) lssormax, tol, wo
     END IF
!!$     ! ... calculate max bandwidth if none is input
!!$     IF (mbw > 0) THEN
!!$        WRITE(fuclog,9185) 'user input', mbw
!!$     ELSE IF (nx >= nz) THEN
!!$        mbw = 2*(2*nz+1) + 1
!!$        WRITE(fuclog,9185) '4*NZ+3', mbw
!!$     ELSE
!!$        mbw = 2*(2*nx+1) + 1
!!$        WRITE(fuclog,9185) '4*NX+3', mbw
!!$9185    FORMAT (10X, 'Matrix maximum bandwidth (', a, '): ',i5)
!!$     END IF
  CASE (2)
     WRITE(fupdef,9235) '--- GMRES Linear Equation Solver selected'
     WRITE(fuclog,9235) '--- GMRES Linear Equation Solver selected'
     READ(fuins,*) m_save_dir, stop_tol_gmres, maxit_gmres, ilu_method, lev_fill, drop_tol
     ! ... set defaults as necessary
     IF(maxit_gmres == 0) maxit_gmres = 100
     IF(m_save_dir == 0) m_save_dir = 5
     IF(stop_tol_gmres == 0._kdp) stop_tol_gmres = 1.e-10_kdp
     IF(ilu_method == 0) ilu_method = 1
     IF(drop_tol == 0._kdp) drop_tol = 0.001_kdp
     WRITE(fupdef,9235) '--- Generalized Minimum-Residual Solver Control ---'
     WRITE(fupdef,9041)  &
          'Number of saved directions ........... ',m_save_dir,  &
          'Convergence tolerance for residual ... ',stop_tol_gmres,  &
          'Maximum number of iterations ......... ',maxit_gmres
9041 FORMAT (tr10,a,i10/tr10,a,1pg10.3/tr10,a,i10)
     SELECT CASE (ilu_method)
     CASE (1)
        IF(lev_fill == 0) lev_fill = 5
        WRITE(fupdef,9141) 'ILU Preconditioner with drop tolerance and fill limit',  &
             'Drop tolerance ....................... ',drop_tol,  &
             'Fill limit, k largest elements ....... ',lev_fill
9141    FORMAT (tr10,a/tr10,a,f10.5/tr10,a,i10)
     CASE (2)
        WRITE(fupdef,9241) 'ILU Preconditioner with level of fill',  &
             'Level of fill ........................ ',lev_fill
9241    FORMAT (tr10,a/tr10,a,i10)
     END SELECT
  END SELECT
  ! ...        ----- weighting & averaging options -----
  ! ...     Ioptupst (logic) - T= upstream weighting; F= midpoint weighting
  ! ...     Potdif - the difference in potential factor: upstream weighting
  ! ...         is applied when the potential difference between two nodes
  ! ...         is greater than this value, otherwise a weighted average is used
  ! ...     Theta - time weighting variable (0 to 1) NOT YET IMPLEMENTED!!!
  IF (informt == 1) READ(fuins,9000)
  READ(fuins,*) ioptupst, potdif, theta
  WRITE(fupdef,9235) '--- Weighting & Averaging Options ---'
  IF (ioptupst) THEN
     WRITE(fupdef,9245)  'Upstream'
  ELSE
     WRITE(fupdef,9245) 'Midpoint'
  END IF
  WRITE(fupdef,9045) potdif, theta
  ! ...        ----- relative permeability data -----
  ! ...     Kodrp:  0 - linear relative permeability
  ! ...             1 - Corey-type relative permeability
  ! ...             2 - facture relative permeability
  ! ...             3 - Corey with Cooley's S(pc) function
  ! ...             4 - van Genuchten/Mualem function
  ! ...     Swr - residual water saturation
  ! ...     Ssr - residual steam saturation
  IF (informt == 1) READ(fuins,9000)
  READ(fuins,*) kodrp, swr, ssr
  ! ... check for errors
  IF (kodrp < 0 .OR. kodrp > 4) THEN
     ! ... Unknown relative permeability function
     WRITE(fustdout,9250) '*** ERROR - Invalid index for relative '//  &
          'permeability function  ***'
     WRITE(fuclog,9250) '*** ERROR - Invalid index for relative '//  &
          'permeability function  ***'
9250 FORMAT (/10X, a)
     ierr(20) = .TRUE.
  END IF
  WRITE(fupdef,9235) '--- Relative Permeability ---'
  SELECT CASE (kodrp)
  CASE (0)
     WRITE(fupdef,9255) 'Linear'
  CASE (1)
     WRITE(fupdef,9255) 'Corey-type'
  CASE (2)
     WRITE(fupdef,9255) 'Fracture'
  CASE (3)
     WRITE(fupdef,9255) 'Corey-type with Cooley Saturation, S(pc)'
  CASE (4)
     WRITE(fupdef,9255) 'van Genuchten-Mualem'
  END SELECT
  WRITE(fupdef,9050) swr, ssr
  ! ...                  ----- rock properties -----
  ! ...   Phfwt - heat capacity (specific heat)
  ! ...   Df - rock density
  ! ...   Beta - rock compressibility
  ! ...   Grav - gravitational constant
  ! ...   Initphi - True - read initial porosity and pressure from file; False - no read
  CALL alloc_parameters(a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: gdata, alloc_parameters"  
     ierr(199) = .TRUE.
     RETURN
  ENDIF
  IF (informt == 1) READ(fuins,9000)
  READ(fuins,*) phfwt(1), df(1), beta(1), grav, initphi
  ! ... open and read initial P and PHI input file if required
  CALL alloc_ic(a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: gdata, alloc_ic"  
     ierr(199) = .TRUE.
     RETURN
  ENDIF
  IF (initphi) THEN
     READ(fuins,9190) phifile
     CALL gfiles(2, phifile, flag)
     READ(fuicpp,9000)
     READ(fuicpp,*) (pinit(ic), ic=1,nxyz)
     READ(fuicpp,9000)
     READ(fuicpp,*) (phiinit(ic), ic=1,nxyz)
  END IF
  ! ...     beginning with version 2.01, specific heat, rock density and
  ! ...       compressibility can vary in space
  ! ...       if uniform - their values may be input here as with earlier
  ! ...       versions, or if a -1 is input their values are input with the
  ! ...       other arrays
  ! ... The value for NARRAY - the number of arrays to input after cell size arrays must be
  ! ...      increased for each rock property as necessary
  narray = 8
  WRITE(fupdef,9235) '--- Rock Properties ---'
  IF (phfwt(1) == -1._KDP) THEN
     WRITE(fupdef,9055) 'Specific Heat of Rock (erg/g-K):   '
     kod(6) = 1
     narray = narray + 1
  ELSE
     DO  ic = 2,nxyz
        phfwt(ic) = phfwt(1)
     END DO
     unitfac(6) = 1._kdp
     unitlabel(6) = '(erg/g-K)'
     WRITE(fupdef,9060) 'Specific Heat of Rock (erg/g-K):  ', phfwt(1)
     kod(6) = 2
  END IF
  IF (df(1) == -1._kdp) THEN
     WRITE(fupdef,9055) 'Rock Density (g/cm^3):             '
     kod(7) = 1
     narray = narray + 1
  ELSE
     DO  ic = 2,nxyz
        df(ic) = df(1)
     END DO
     unitfac(7) = 1._kdp
     unitlabel(7) = '(g/cm^3)'
     WRITE(fupdef,9060) 'Rock Density (g/cm^3):            ', df(1)
     kod(7) = 2
  END IF
  IF (beta(1) == -1._kdp) THEN
     WRITE(fupdef,9055) 'Rock Compressibility (cm^2/dyne):  '
     kod(8) = 1
     narray = narray + 1
  ELSE
     DO  ic = 2,nxyz
        beta(ic) = beta(1)
     END DO
     unitfac(8) = 1._kdp
     unitlabel(8) = '(cm^2/dyne)'
     WRITE(fupdef,9060) 'Rock Compressibility (cm^2/dyne): ',beta(1)
     kod(8) = 2
  END IF
  WRITE(fupdef,9065) 'Gravitational Constant (cm/sec^2): ',grav
  ! ... Error:Original P & PHI not input but should be because BETA > 0
  IF (timestr > 0._kdp .AND. beta(1) > 0._kdp) THEN
     WRITE(fustdout,9260)
     WRITE(fuclog,9260)
9260 FORMAT (/10X, '***** WARNING  *****'/10X,  &
          'Initial Condition P & PHI not provided but should be '//  &
          'because BETA > 0')
  END IF
  IF (initphi) THEN
     WRITE(fupdef,9180) 'Input file for original P & PHI = ',  &
          TRIM(phifile)
     IF (timestr == 0._kdp .OR. beta(1) == 0._kdp) THEN
        WRITE(fustdout,9200) TRIM(phifile)
        WRITE(fuclog,9200) TRIM(phifile)
9200    FORMAT (10X, '***** ERROR *****'/10X,  &
             'Initial condition P & PHI should not be input from file '/  &
             tr20,a/10X, 'because either BETA = 0 or start time = 0')
        ierr(10) = .TRUE.
        !..          return
     END IF
  END IF
  CALL alloc_variables(a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: gdata, alloc_variables"  
     ierr(199) = .TRUE.
     RETURN
  ENDIF
  ! ...  cylindrical coordinates
  ! ...    Irad: True, cylindrical coords.; False, Cartesian coords.
  ! ...    Rsq(1): well radius
  ! ...    Rsq(Nxx): outer radius of region
  IF (informt == 1) READ(fuins,9000)
  READ(fuins,*) irad, rsq(1), rsq(nxx)
  IF (irad) THEN
     WRITE(fupdef,9070) rsq(1), rsq(nxx)
     ny = 1
     IF (informt >= 3) narray = narray - 1        ! no y perm read for cylindrical when informt >= 3
  END IF
  ! ... print/plotting controls
  CALL printopt(99)
  ! ... At first read of print options, open output files for printdata and plotfile
  ufile = ' '
  CALL gfiles(3, ufile, flag)
  ! ...         ----- Node information for each mesh Slice -----
  ! ...    read in a sequence of records for each slice
  ! ...    Ifmt - input format to describe mesh region
  ! ...           0,F = old style format
  ! ...           1,T = version 1.02 style input format
  ! ...           2   = version 2.01 and version 3 style with ROCKTYPES
  ! ...    Nb(J) - number of active blocks per slice
  ! ...    Nout - number of blocks outside problem domain for this slice
  ! ...    Nconst - number of constant P & H nodes in this slice
  WRITE(fupdef,9235) '--- Node Sequence Numbers for the Mesh ---'
  ! ... assign formats for given printing width
  IF (ioptpr(2) == 0) THEN
     ! ...       80 columns
     cfmtc = cfmtc1
     cfmtd = cfmtd1
  ELSE IF (ioptpr(2) == 1) THEN
     ! ...       132 columns
     cfmtc = cfmtc2
     cfmtd = cfmtd2
  ELSE IF (ioptpr(2) >= 2) THEN
     ! ...       159 columns
     cfmtc = cfmtc3
     cfmtd = cfmtd3
  END IF
  CALL alloc_bc(a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: gdata, alloc_bc"  
     ierr(199) = .TRUE.
     RETURN
  ENDIF
  ! ... Input node definition information by x-z slice
  nodseqnoglb = 0          ! ... initialize global node number
  DO  j=1,ny
     ! ...        for old input format FORM2
     SELECT CASE (informt)
     CASE (1)
        READ(fuins,9000)
        READ(fuins,*) newfmt, nb(j), nout, nconst
        ifmt = 0
        IF (newfmt) ifmt = 1
     CASE (2)
        READ(fuins,*,ERR=380) ifmt, nb(j), nout, nconst
        ! ... for input format VER2, if a T or F is read in place for Ifmt
        ! ... (as was required for Vers 1.02), an i/o error will occur,
        !     and execution to go to statment label 380 to write an error message
     CASE (3)
        ifmt = 2
     END SELECT
     ! ...     for old and new formats (Ifmt=0,1)
     ! ...          check to see that Nb(j) + Nout + Nconst sum to nxz
     ! ...     this is not necessary for latest (ROCKTYPE) Ifmt=2 format
     IF (ifmt < 2 .AND. nb(j)+nout+nconst /= nxz) THEN
        ! ... Problem with number of nodes in slice
        WRITE(fustdout,9265) j
        WRITE(fuclog,9265) j
9265    FORMAT ('ERROR in number of nodes in slice',i5)
        ierr(35) = .TRUE.
     END IF
     ! ...         ######### old input format ########
     ! ...           for locations of active, boundary, and outside nodes
     ! ...               read sequencing of blocks in slice
     ! ...            Np(ik,j) - sequence numbering of block I,K  / slice
     ! ...              NOTE:  Read zeros for missing nodes (not part of domain)
     ! ...              NOTE:  Read -1 for constant pressure-enthalpy nodes
     IF (ifmt == 0) THEN
        READ(fuins,*) (np(ik,j), ik=1,nxz)
     ELSEIF (ifmt == 1) THEN
        ! ...         ######### Version 1.02 format ########
        ! ...             for locations of active, boundary, and outside nodes
        ! ...             for array Np(ik,j)
        ! ...             assign -1 to constant nodes and
        ! ...             temporarily assign -99 to nodes outside of the domain
        ! ...                   then later reassign a 0 to these nodes
        DO  ldum = 1,nconst
           READ(fuins,*) i, k, lll
           ik = (k-1)*nx + i
           np(ik,j) = -1           ! ... constant P & H nodes
           ibc(ik,j) = 11
           ! ...            nodes outside domain to the right of the constant P & H node
           IF (lll > 0 .AND. i < nx) THEN
              DO  ii = i + 1,nx
                 ik = (k-1)*nx + ii
                 np(ik,j) = -99
                 ibc(ik,j) = -1
              END DO
           END IF
           ! ...            nodes outside domain to the left of the constant P & H node
           IF (lll < 0 .AND. i > 1) THEN
              DO  ii = 1,i - 1
                 ik = (k-1)*nx + ii
                 np(ik,j) = -99
                 ibc(ik,j) = -1
              END DO
           END IF
        END DO
     ELSEIF (ifmt == 2) THEN
        ! ...         ######### Version 2.1 and 3 ROCKTYPE format ########
        ! ... Identification of active nodes with rock type, specified p & h boundary
        ! ...      nodes, seepage boundary nodes, and excluded nodes.
        ! ... read TOP or BOTTOM to determine whether rows are
        ! ... input from the top first or bottom first
        READ(fuins,'(A)') ustring
        ustring = ADJUSTL(ustring)
        ustring = uppercase(ustring(1:10))
        IF (ustring(1:3) /= 'TOP' .AND. ustring(1:3) /= 'BOT') THEN
           WRITE(fustdout,9235) 'SLICE Input requires TOP or BOTTOM keyword'
           WRITE(fuclog,9235) 'SLICE Input requires TOP or BOTTOM keyword'
9235       FORMAT (/a)
           ierr(36) = .TRUE.
        END IF
        ! ... read rocktype indices
        ! ...              0 = out of problem domain
        ! ...             -n = specified pressure and enthalpy node of rock type n
        ! ...              n = active node of rock type n
        ! ...            -1n = specified pressure node with associated enthalpy or temperature
        ! ...                  of rock type n
        ! ...            -3n = seepage node of rock type n
        ! ...            -5n = combination seepage node and precipitation recharge node 
        ! ...                  of rock type n
        ! ...  Rock type is a 1 or 3 digit integer depending on the input format flag
        mprtp1 = 10
        IF(informtp >= 1) mprtp1 = 1000
        IF (ustring(1:3) == 'TOP') THEN
           DO  k = nz,1,-1
              icxf = (k-1)*nx + 1 + (j-1)*nxz
              icxl = icxf + nx - 1
              READ(fuins,*) (icrxtype(ic), ic=icxf,icxl)
           END DO
        ELSE
           icxf = 1 + (j-1)*nxz
           icxl = j*nxz
           READ(fuins,*) (icrxtype(ic), ic=icxf,icxl)
        END IF
        ! ... determine the number of active nodes and constant nodes
        ! ... and for array Np(ik,j) assign -1 to constant nodes and
        ! ... temporarily assign -99 to nodes outside of the domain
        ! ... then later reassign a 0 to these nodes
        ! ... Seepage cells are assigned the rock type index
        ! ... IBC array stores boundary condition information
        nconst = 0
        nspecpbc = 0
        nseep = 0
        nout = 0
        nb(j) = 0
        DO  k = 1,nz
           DO  i = 1,nx
              ik = (k-1)*nx + i
              ic = (k-1)*nx + i + (j-1)*nxz
              irx = MOD(ABS(icrxtype(ic)),mprtp1)
              ibctx = icrxtype(ic)/mprtp1
              IF (icrxtype(ic) > 0) THEN
                 ! ... active nodes    - temporarily set np=0
                 np(ik,j) = 0
                 ibc(ik,j) = 0
                 nb(j) = nb(j) + 1
              ELSE IF (icrxtype(ic) < 0 .and. ibctx == 0) THEN
                 ! ... specified p and h value nodes 
                 np(ik,j) = -1
                 ibc(ik,j) = 11
                 nconst = nconst + 1
              ELSE IF (ibctx == -1) THEN
                 ! ... specified p with associated h or t nodes
                 np(ik,j) = 0
                 ibc(ik,j) = 10
                 nspecpbc = nspecpbc + 1
                 nb(j) = nb(j) + 1
              ELSE IF (ibctx == -3) THEN
                 ! ... seepage nodes
                 np(ik,j) = 0
                 ibc(ik,j) = 30
                 nseep = nseep + 1
                 nb(j) = nb(j) + 1
              ELSE IF (ibctx == -5) THEN
                 ! ... combination seepage and precipitation nodes
                 np(ik,j) = 0
                 ibc(ik,j) = 50
                 nseep = nseep + 1
                 nb(j) = nb(j) + 1
              ELSE
                 ! ... out of domain nodes  - temporarily set np=-99
                 np(ik,j) = -99
                 ibc(ik,j) = -1
                 nout = nout + 1
              END IF
              icrxtype(ic) = irx
              !   From here onward, ircxtype is non-negative and less than mprtp1 (10 or 1000, depending on version). 
              !   That is, icrxtype is the rock type index with no BC information. BC info stored in array ibc.
           END DO
        END DO
     END IF
     IF (ifmt == 1 .OR. ifmt == 2) THEN
        ! ... For SSOR solver, for array np(ik,j) assign sequence numbers to active nodes
        ! ...       (not specified p,h nodes) and assign 0 to nodes outside domain
        ! ... For GMRES iterative solver, assign node sequence numbers to active nodes over
        ! ...      entire mesh. Use natural order of sweeping the mesh.
        nodseqno = 0             ! ... initialize node sequence number
        seeping = .FALSE.        ! ... initialize seepage flag
        IF(slmeth == 1) THEN
           IF (nz >= nx) THEN    ! ...  number the nodes row by row
              DO  k = 1,nz
                 DO  i = 1,nx
                    ik = (k-1)*nx + i
                    ic = (k-1)*nx + i + (j-1)*nxz
                    IF (np(ik,j) >= 0) THEN
                       nodseqno = nodseqno + 1
                       np(ik,j) = nodseqno
                       !$$                    nodseqnoglb = nodseqnoglb + 1
                       !$$                    mrno(ic) = nodseqnoglb
                    ELSE
                       IF (np(ik,j) == -99) np(ik,j) = 0
                    END IF
                 END DO
              END DO
           ELSE        ! ... number the nodes column by column
              DO  i = 1,nx
                 DO  k = 1,nz
                    ik = (k-1)*nx + i
                    IF (np(ik,j) >= 0) THEN
                       nodseqno = nodseqno + 1
                       np(ik,j) = nodseqno
                    ELSE
                       IF (np(ik,j) == -99) np(ik,j) = 0
                    END IF
                 END DO
              END DO
           END IF
        ELSEIF(slmeth ==2) THEN     ! ...  number the nodes row by row
           DO  k = 1,nz
              DO  i = 1,nx
                 ik = (k-1)*nx + i
                 ic = (k-1)*nx + i + (j-1)*nxz
                 IF (np(ik,j) /= -99) THEN
                    nodseqno = nodseqno + 1
                    np(ik,j) = nodseqno
                    nodseqnoglb = nodseqnoglb + 1
                    mrno(ic) = nodseqnoglb
                 ELSE
                    np(ik,j) = 0
                 END IF
              END DO
           END DO
        END IF
     END IF
     ! ... write node numbers for either input format
     IF(slmeth ==1) THEN
        WRITE(fupdef,9075) j, nb(j), nconst, nout
   9075 FORMAT (/tr7,'Node Numbers for SLICE',i4/8X,i4,  &
          ' Active Nodes;',i4, ' Constant Value Nodes (-1);',i4,  &
          ' Inactive Nodes (0)')
        WRITE(fupdef,cfmtc) (i, i=1,nx)
        WRITE(fupdef,9190) ' '
        DO  k = nz,1,-1
           ikxf = (k-1)*nx + 1
           ikxl = k*nx
           WRITE(fupdef,cfmtd) k, (np(ik,j), ik=ikxf,ikxl)
        END DO
     ELSEIF(slmeth ==2) THEN
        WRITE(fupdef,9074) j, nb(j), nconst, nout
   9074 FORMAT (/tr7,'Node Numbers for SLICE',i4/8X,i4,  &
          ' Active Nodes;',i4, ' Constant Value Nodes;',i4,  &
          ' Inactive Nodes (0)')
        WRITE(fupdef,cfmtc) (i, i=1,nx)
        WRITE(fupdef,9190) ' '
        DO  k = nz,1,-1
           icxf = (k-1)*nx + 1 + (j-1)*nxz
           icxl = icxf + nx - 1
           WRITE(fupdef,cfmtd) k, (mrno(ic), ic=icxf,icxl)
        END DO
     ENDIF
     ! ... write b.c. index numbers for either input format
     WRITE(fupdef,9076) 'Boundary Condition Index Numbers for SLICE',j,  &
          nconst,' Constant Value P & H & T Nodes (11);',  &
          nspecpbc,' Specified P, associated H and T Nodes (10);',  &
          nseep,' Seepage Surface Nodes (30 or 50);',  &
          nout,' Inactive External Nodes (-1)'
9076 FORMAT (/tr7,a,i4/tr8,i4,a,i4,a/tr8,i4,a,i4,a)
     WRITE(fupdef,cfmtc) (i, i=1,nx)
     WRITE(fupdef,9190) ' '
     DO  k = nz,1,-1
        ikxf = (k-1)*nx + 1
        ikxl = k*nx
        WRITE(fupdef,cfmtd) k, (ibc(ik,j), ik=ikxf,ikxl)
     END DO
  END DO     ! ... end of loop for slice J
  ! ... Assign the i,j,k, indices to array npp for only the active and seepage
  ! ...       nodes
  ! ... and assign ic indices to array npic
  ! ... array elements 1 - npicmax are indices of active and seepage nodes and
  ! ... array elements npicmax+1 - npiccons are indices of constant P and H value nodes
  lll = 0
  mrnomax = 0
  DO  j=1,ny
     DO  k=1,nz
        DO  i=1,nx
           ik = (k-1)*nx + i
           ! ... nn and np(ik,j) are the sequence numbers of cells
           ! ... ml are the sequence numbers of the active and seepage cells for a given
           ! ...      slice
           ! ...  -1:specified p,h node, 0:node out of domain
           ml = np(ik,j)
           IF (ml > 0) THEN
              ! ... Store pointers for SSOR linear equation solver array and count 
              ! ...      active and seepage cells
              npp(2*ml-1, j) = i
              npp(2*ml, j) = k
              lll = lll + 1
              ic = (k-1)*nx + i + (j-1)*nxz
              npic(lll) = ic
              mrnomax = MAX(mrnomax,mrno(ic))
           END IF
        END DO
     END DO
  END DO
  npicmax = lll
  DO  j=1,ny
     DO  k=1,nz
        DO  i=1,nx
           ik = (k-1)*nx + i
           IF (np(ik,j) <= -1) THEN
              ! ... Store constant value nodes in npic and count
              ic = (k-1)*nx + i + (j-1)*nxz
              lll = lll + 1
              npic(lll) = ic
           END IF
        END DO
     END DO
  END DO
  npiccons = lll
  ! ... Check for errors in node sequence numbers
  DO  j=1,ny
     nby = nb(j)
     DO  ml=1,nby-1
        mlnext = ml + 1
        ! ... determine the I and K indices of adjacent active nodes
        i = npp(2*ml-1,j)
        k = npp(2*ml,j)
        ik = (k-1)*nx + i
        inext = npp(2*mlnext-1,j)
        knext = npp(2*mlnext,j)
        iknext = (knext-1)*nx + inext
        IF (np(ik,j)+1 /= np(iknext,j)) THEN
           WRITE(fustdout,9205) np(ik,j)
           WRITE(fuclog,9205) np(ik,j)
9205       FORMAT (/'***  ERROR in node sequence numbers ***'/  &
                TR10, 'Check nodes near sequence number',i5)
           ierr(1) = .TRUE.
        END IF
     END DO
  END DO
  WRITE(fuclog,9235) '--- Solver Information ---'
  !   ------- Check to see if matrix bandwidth "Mbw" is large enough -------
  ! ...   for way the nodes are ordered, when the old input format is used
  IF (ifmt == 0 .AND. nz > 1 .AND. nx > 1) THEN
     ! ... determine whether nodes are numbered by rows or columns
     lrowcol = 0
     lll = 0
     ndumy = 0
     DO
        lll = lll + 1
        ic = npic(lll)
        j = IDINT(DBLE(ic)/DBLE(nxz)) + 1
        i = npp(2*lll-1, j)
        k = npp(2*lll, j)
        ik = (k-1)*nx + i
        IF (ic >= nxyz) ndumy = 1
        ! ... check if adjacent nodes are active nodes
        IF (nb(j) <= 1) THEN
           ! ... if there is only 0 or 1 active nodes
           lrowcol = 3
           ndumy = 1
        ELSE IF (i /= nx .AND. k /= nz .AND. np(ik+1,j) > 0 .AND.  &
             np(ik+nx,j) > 0) THEN
           ! ... check for row by row numbering
           IF (k == npp(2*(lll+1),j)) lrowcol = 1
           ! ... check for column by column numbering
           IF (i == npp(2*(lll+1)-1,j)) lrowcol = 2
           ndumy = 1
        END IF
        IF (ndumy == 1) EXIT
     END DO
     ! ... check if row or column determination was successful
     ! ... and write warning about node numbering if necessary
     IF (lrowcol == 0) THEN
        WRITE(fustdout,'(/a/tr5,a/tr5,a)')   &
             '*** WARNING - Problem with Node Checking ***',  &
             'Unable to determine whether nodes are numbered by row '//  &
             'or column',  &
             'Check that node numbering and matrix bandwidth are correct'
        WRITE(fuclog,'(/a/tr5,a/tr5,a)')   &
             '*** WARNING - Problem with Node Checking ***',  &
             'Unable to determine whether nodes are numbered by row '//  &
             'or column',  &
             'Check that node numbering and matrix bandwidth are correct'
     ELSEIF (lrowcol == 1) THEN
        mbw = 2*(2*nx+1) + 1
        WRITE(fuclog,9185) '4*NX+3', mbw
9185    FORMAT (10X, 'Matrix maximum bandwidth (', a, '): ',i5)
     ELSEIF (lrowcol == 2) THEN
        mbw = 2*(2*nz+1) + 1
        WRITE(fuclog,9185) '4*NZ+3', mbw
     ENDIF
  END IF
  IF(slmeth == 1) THEN
     IF (nx >= nz) THEN
        mbw = 2*(2*nz+1) + 1
        WRITE(fuclog,9185) '4*NZ+3', mbw
     ELSE
        mbw = 2*(2*nx+1) + 1
        WRITE(fuclog,9185) '4*NX+3', mbw
     END IF
     CALL alloc_solver(a_err)
     IF (a_err /= 0) THEN  
        PRINT *, "Array allocation failed: gdata, alloc_solver"
        ierr(199) = .TRUE.
        RETURN
     ENDIF
  ELSEIF(slmeth == 2) THEN
     CALL ldci       ! ... Set connection list
     CALL ldiaja     ! ... Determine size and structure of A array in CSR storage
     CALL alloc_solver_gmres(a_err)
     IF (a_err /= 0) THEN  
        PRINT *, "Array allocation failed: gdata, alloc_solver_gmres"
        ierr(199) = .TRUE.
        RETURN
     ENDIF
     neq = 2*mrnomax
     SELECT CASE (ilu_method)
     CASE (1)          ! ... drop tolerance with fill limit
        nnz_ilu = neq*(2*lev_fill+1)
        CALL alloc_ilu_gmres(a_err)
        WRITE(fuclog,9001) 'Number of non-zero elements in A ........... ', ia(neq+1) - ia(1)
        WRITE(fuclog,9001) 'Number of non-zero elements in ILU(k,tol) .. ', nnz_ilu
9001    FORMAT(tr10,a,i8)
     CASE (2)          ! ... level k fill
        CALL fill_stor(neq,ja,ia,lev_fill,nnz_ilu)
        CALL alloc_ilu_gmres(a_err)
        WRITE(fuclog,9001) 'Number of non-zero elements in A ........... ', ia(neq+1) - ia(1)
        WRITE(fuclog,9001) 'Number of non-zero elements in ILU(k) ...... ', nnz_ilu
     END SELECT
     IF (a_err /= 0) THEN  
        PRINT *, "Array allocation failed: gdata, alloc_ilupc"
        ierr(199) = .TRUE.
        RETURN
     ENDIF
  ENDIF
  IF(informt == 3) THEN
     ! ... Version 3 data
     ! ... Read flag for unconfined flow, atmospheric pressure
     READ(fuins,*) unconfined,patm
     ! ...      saturation-pressure parameter
     IF(unconfined) THEN
        WRITE(fupdef,9902)  &
             'Unconfined flow is allowed in the region for this simulation',  &
             'Atmospheric pressure (dyne/cm^2) .......',patm
9902    FORMAT(/tr5,a/2(tr5,a,1PG14.5))
        IF(kodrp == 0) THEN
           READ(fuins,*) pwb,pwr
           WRITE(fupdef,9903)  &
                'Bubble point pressure (Linear parameter) (dyne/cm^2) ...', pwb,  &
                'Residual pressure (Linear parameter) (dyne/cm^2) .......', pwr
        ELSE IF(kodrp == 1) THEN
           READ(fuins,*) pwb
           ! ***... add lambda to this
           WRITE(fupdef,9903)  &
                'Bubble point pressure (Corey parameter) (dyne/cm^2) ....', pwb
        ELSE IF(kodrp == 3) THEN
           READ(fuins,*) pwb,bb,cc
           WRITE(fupdef,9903)  &
                'Bubble point pressure (Cooley parameter) (dyne/cm^2) ....', pwb,  &
                'Brutsaert parameter C (Cooley parameter A) '//  &
                '(dyne/cm^2)^cc ....', bb,  &
                'Brutsaert exponent parameter c (Cooley parameter c) '// '(-) ....',cc
!!$        ELSE IF(kodrp == 4) THEN
!!$           READ(fuins,*) pwb, lambda
!!$           gamma = 1._kdp -1._kdp/lambda 
!!$           WRITE(fupdef,9903)  &
!!$                'Scaling pressure (van Genuchten parameter) (dyne/cm^2) ....', pwb, &
!!$                'Pore size index (van Genuchten parameter) (-) .............', lambda, &
!!$                'Inflection point index (calculated van Genuchten parameter) (-) ......', gamma
        END IF
9903    FORMAT(/3(tr5,a,1PG14.5/))
        ! ... Increase number of arrays to be read to include precipitation
        narray = narray + 2
     END IF
  END IF
  ! ...  read Dx,Dy,Dz arrays first so the cell sizes
  ! ...     can be used to calculate the distribution of
  ! ...     other spatial properties
  IF (informt == 1) THEN
     READ(fuins,9000)
     READ(fuins,9000)
  END IF
  IF(informt >= 3 .AND. irad) THEN     ! if informt <= 2, still have to read Y spacing
     WRITE(fuclog,9235) '--- Read R and Z dimensions of cells ---'
     CALL readarry(2)
  ELSE
     WRITE(fuclog,9235) '--- Read X, Y, and Z dimensions of cells ---'
     CALL readarry(3)
  END IF
  ! ... assign format for writing node and cell spacing arrays depending on screen display width
  IF (ioptpr(2) == 0) THEN
     cfmta = cfmta1
     cfmtb = cfmtb1
  ELSE IF (ioptpr(2) == 1) THEN
     cfmta = cfmta2
     cfmtb = cfmtb2
  ELSE IF (ioptpr(2) >= 2) THEN
     cfmta = cfmta3
     cfmtb = cfmtb3
  END IF
  ! ...     -------- write input and calculated values for arrays --------
  IF (irad) THEN
     ! ... XSPACING of nodes for cylindrical coordinates
     IF (dx(1) == -1._kdp) THEN          ! ...      automatic node spacing and cell sizes
        WRITE(fupdef,9275) '--- Input radii for Cylindrical Coordinates ---',  &
             'Radii of nodes automatically computed using input:',  &
             'Radius of first node:                 ',xx(1)/unitfac(12),  &
             TRIM(unitlabel(12)),  & 
             'Maximum radius of simulation region:  ',rsq(nxx)/unitfac(12),  &
             TRIM(unitlabel(12)),  & 
             'Number of nodes in R directon:        ',nx
9275    FORMAT (/a/tr10,a/tr15,a,1pg10.3,tr1,a/tr15,a,1pg10.3,tr1,a/tr15,a,i6)
        WRITE(fupdef,9235)  &
             'Node locations in the R-direction '//TRIM(unitlabel(12))
        WRITE(fupdef,cfmtb) 'Column Numbers (I)',(i, i=1,nx)
        WRITE(fupdef,cfmta) (xx(i)/unitfac(12), i=1,nx)
        ! ... calculate and write squared boundaries of cells
        rsq(1) = rsq(1)**2
        rsq(nxx) = rsq(nxx)**2
        DO  i = 2,nx
           rsq(i) = (xx(i)**2-xx(i-1)**2)/LOG(xx(i)**2/xx(i-1)**2)
        END DO
        WRITE(fupdef,9280) 'Cylindrical Coordinates ---',  &
             'Note: grid is discretized in r^2 for cylindrical coordinates',  &
             'First cell bounded by well radius squared:  ',rsq(1)/(unitfac(12)*unitfac(12)),  &
             TRIM(unitlabel(12))//'^2',  &
             'Input value for well radius:                ',SQRT(rsq(1))/unitfac(12),  &
             TRIM(unitlabel(12))
9280    FORMAT (/a/tr10,a/tr15,a,1pg10.3,tr1,a/tr15,a,1pg10.3,tr1,a)
        WRITE(fupdef,'(/a)') '--- Cell boundary R-locations squared '//  &
             TRIM(unitlabel(12))//'^2 ---'
        WRITE(fupdef,cfmta) (rsq(i)/(unitfac(12)*unitfac(12)), i=1,nxx)
        ! ... calculate and write sizes of cells
        DO  i = 1,nx
           dx(i) = SQRT(rsq(i+1)) - SQRT(rsq(i))
        END DO
        WRITE(fupdef,'(/a/tr10,a)')  &
             '--- R-dimensions of cylindrical cells '//TRIM(unitlabel(12))//' ---',  &
             'Calculated from cell boundary locations squared'
        WRITE(fupdef,cfmtb) 'Column Numbers (I)', (i, i=1,nx)
        WRITE(fupdef,cfmta) (dx(i)/unitfac(12), i=1,nx)
     ELSE                ! ...  cell radial sizes are individually specified
        WRITE(fupdef,'(/a/tr10,a/tr10,a)')  &
             '--- Input values for R-dimensions of cylindrical cells '//  &
             TRIM(unitlabel(12))//' ---',  &
             '***Cell size is recalculated from cell boundary locations squared',  &
             '   ***   Specified maximum aquifer radius not used'
        WRITE(fupdef,cfmtb) 'Column Numbers (I)', (i, i=1,nx)
        WRITE(fupdef,cfmta) (dx(i)/unitfac(12), i=1,nx)
        ! ... calculate locations of nodes and squared boundaries of cells
        ! ... node locations are at radial midpoint of cell
        xx(1) = rsq(1) + 0.5_kdp*dx(1)  
        DO  i = 2,nx
           xx(i) = xx(i-1) + 0.5_kdp*(dx(i-1) + dx(i))
           rsq(i) = rsq(i-1) + dx(i-1)          ! ... cell face radii at present
        END DO
        WRITE(fupdef,9285) 'Node locations in R-direction '//TRIM(unitlabel(12)),  &
             'Calculated from input cell sizes and well radius'
9285    FORMAT (/tr5,a/tr10,a)
        WRITE(fupdef,cfmta) (xx(i)/unitfac(12), i=1,nx)
        DO  i = 1,nxx
           rsq(i) = rsq(i)**2
        END DO
        WRITE(fupdef,'(/tr5,a/tr10,a)') '--- Cell boundary R-locations squared '//  &
             TRIM(unitlabel(12))//'^2 ---'
        WRITE(fupdef,cfmta) (rsq(i)/(unitfac(12)*unitfac(12)), i=1,nxx)
     END IF
     ! ... Print out the distance along R axis to sides of cells in user units
     WRITE(fupdef,9195) 'Cell face locations in R-direction '//TRIM(unitlabel(12))
     WRITE(fupdef,cfmta) (SQRT(rsq(i))/unitfac(12), i=1,nxx)
  ELSE                ! ... XSPACING of nodes for Cartesian coords
     WRITE(fupdef,9235)  &
          '--- X-dimensions of cells for Cartesian Coordinates '//TRIM(unitlabel(12))//' ---'
     IF (kod(12) == 2) THEN          ! ... uniform X dimensions of cells
        WRITE(fupdef,cfmtb) 'Column Numbers (I)', (i, i=1,nx)
        WRITE(fupdef,9130) dx(1)/unitfac(12)
     ELSE                            ! ... variable X dimensions of cells
        WRITE(fupdef,cfmtb) 'Column Numbers (I)', (i, i=1,nx)
        WRITE(fupdef,cfmta) (dx(i)/unitfac(12), i=1,nx)
     END IF
     ! ... compute distance along X to node location, XX
     ! ...     (left hand side of grid is zero datum)
     xx(1) = 0.5_kdp*dx(1)
     DO  i = 2,nx
        xx(i) = xx(i-1) + 0.5_kdp*(dx(i-1) + dx(i))
     END DO
     ! ... write distance along X axis to cell node location in user units
     WRITE(fupdef,9195) 'Node locations in X-direction '//TRIM(unitlabel(12))
     WRITE(fupdef,cfmta) (xx(i)/unitfac(12), i=1,nx)
     ! ... Print out the distance along X axis to sides of cells in user units
     WRITE(fupdef,9195) 'Cell face locations in X-direction '//TRIM(unitlabel(12))
     WRITE(fupdef,cfmta) (xx(1)-0.5_kdp*dx(1))/unitfac(12),  &
          ((xx(i)+0.5_kdp*dx(i))/unitfac(12), i=1,nx)
  END IF
  ! ...                          YSPACING of nodes
  IF (irad) THEN
     dy(1) = 1._kdp
     unitfac(13) = unitfac(12)         ! ... make y direction same units as r direction
     unitlabel(13) = unitlabel(12)     ! ...     only used for generic 3d plot output
  ELSE
     WRITE(fupdef,9235) '--- Y-dimensions of cells '//TRIM(unitlabel(13))//' ---'
     ! ... Kod(13)=1 for variable Y dimensions of cells; =2 for constant Y dimensions
     IF (kod(13) == 2) THEN
        WRITE(fupdef,9130) dy(1)/unitfac(13)
     ELSE
        WRITE(fupdef,cfmta) (dy(j)/unitfac(13), j=1,ny)
     END IF
  END IF
  ! ... compute distance along Y to node location, YY
  ! ...   (front cell boundary plane of grid is zero datum)
  yy(1) = 0.5_kdp*dy(1)
  DO  j = 2,ny
     yy(j) = yy(j-1) + 0.5_kdp*(dy(j-1) + dy(j))
  END DO
  IF (ny > 1) THEN
     WRITE(fupdef,9195) 'Node locations in Y-direction '//TRIM(unitlabel(13))
     WRITE(fupdef,cfmta) (yy(j)/unitfac(13), j=1,ny)
     ! ... Print out the distance along Y axis to rear side of cells in user units
     WRITE(fupdef,9195) 'Cell face locations in Y-direction '//TRIM(unitlabel(13))
     WRITE(fupdef,cfmta) (yy(1)-0.5_kdp*dy(1))/unitfac(13),  &
          ((yy(j)+0.5_kdp*dy(j))/unitfac(13), j=1,ny)
  END IF
  ! ...                     ZSPACING of nodes
  WRITE(fupdef,9235) '--- Z-dimensions of cells ---'
  ! ... Kod(14)=1 for variable Z dimensions of cells; =2 for constant Z dimensions
  ! ...      Compute Zz - elevation of the node point at center of a cell
  ! ...      Compute Zzz - the distance between a node and the one above it
  ! ...         (lowest face of grid is zero datum)
  zz(1) = 0.5_kdp*dz(1)
  DO  k = 2,nz
     zz(k) = zz(k-1) + 0.5_kdp*(dz(k-1) + dz(k))
     zzz(k-1) = zz(k) - zz(k-1)
  END DO
  WRITE(fupdef,'(tr22,A,A/tr11,A,A/tr14,A,tr13,A,tr13,a,tr19,a,tr14,a)')  &
       'Elevation (datum = region cell bottom)   ',  &
       'Depth (datum = region cell top)',  &
       'Cell Size      node location    cell top             ',  &
       'node location    cell bottom',  &
       TRIM(unitlabel(14)), TRIM(unitlabel(14)), TRIM(unitlabel(14)),  &
       TRIM(unitlabel(14)), TRIM(unitlabel(14))
  WRITE(fupdef,9085) (k, dz(k)/unitfac(14),zz(k)/unitfac(14),  &
       (zz(k)+0.5_kdp*dz(k))/unitfac(14),  &
       ((zz(nz)+0.5_kdp*dz(nz))-zz(k))/unitfac(14),  &
       ((zz(nz)+0.5_kdp*dz(nz))-(zz(k)-0.5_kdp*dz(k)))/unitfac(14),  &
       k=nz,1,-1)
9085 FORMAT ('Row',i4,1pe14.4,1pe16.4,1pe15.4,tr8,1pe15.4,1pe17.4)
  ! ... Determine the K value and the elevation of the top node in
  ! ...     each column that is in the active region
  ! ...     Also calculate and store the z-face area
  DO  j = 1,ny
     loop320:  DO  i = 1,nx
        ij = (j-1)*nx + i
        DO  k = nz,1,-1
           ik = (k-1)*nx + i
           IF (np(ik,j) /= 0) THEN
              ktop(ij) = k
              zztop(ij) = zz(k)
              zls(ij) = zz(k) + 0.5_kdp*dz(k)
              land_surf_area(ij) = dx(i)*dy(j)
              CYCLE loop320
           END IF
        END DO
     END DO loop320
  END DO
  ! ... determine average node spacing in X, Y, and Z
  avgdx = (xx(nx)+0.5_kdp*dx(nx))/REAL(nx,KIND=kdp)
  avgdy = (yy(ny)+0.5_kdp*dy(ny))/REAL(ny,KIND=kdp)
  avgdz = (zz(nz)+0.5_kdp*dz(nz))/REAL(nz,KIND=kdp)
  ! ... determine which ROCK types are valid and the number of types
  IF (ifmt == 2) THEN
     DO  irx=1,mprock
        ilrxtype(irx) = .FALSE.
     END DO
     DO  icdum = 1,npiccons
        ic = npic(icdum)
        ilrxtype(icrxtype(ic)) = .TRUE.
     END DO
  END IF
  nrxtype = 0
  DO  irx = 1, mprock
     IF (ilrxtype(irx)) THEN
        nrxtype = nrxtype + 1
        irxused(nrxtype) = irx
     END IF
  END DO
  ! ...      ----- read Phi,Xk,Zk,Yk,Xkc,P, and H or T, arrays -----
  ! ... initialize values for RXPARM
  ! ...   if values for specific heat, rx density, or compressibility
  ! ...   were input above then assign those values here
  DO  lll = 1,nrxtype
     irx = irxused(lll)
     DO  kodparm = 1, mprxparm
        rxparm(kodparm,irx) = -99.999_KDP
        IF (kodparm == 6 .AND. kod(6) == 2) rxparm(kodparm,irx) = phfwt(1)
        IF (kodparm == 7 .AND. kod(7) == 2) rxparm(kodparm,irx) = df(1)
        IF (kodparm == 8 .AND. kod(8) == 2) rxparm(kodparm,irx) = beta(1)
     END DO
  END DO
  WRITE(fupdef, '(/A,I2,A)') '--- Read the remaining ', narray, ' arrays ---'
  WRITE(fuclog, '(/A,I2,A)') '--- Read the remaining ', narray, ' arrays ---'
  CALL readarry(narray)
  DO ie=1,200
     IF(ierr(ie)) THEN
        RETURN          ! ... abort if array reads are in error
     END IF
  END DO
  ! ... PRESSURE and ENTHALPY or TEMPERATURE 
  CALL inputph(.TRUE.)
  ! ... Basal heat flux (conduction) and precipitation flux
  CALL init_wr_bc(.TRUE.)
  ! ... write initial parameters for rock types
  ! ...  POROSITY, PERMEABILITY, THERMAL CONDUCTIVITY SPECIFIC HEAT,
  ! ...       ROCK DENSITY, COMPRESSIBILITY
  CALL init_wr_pmp(.TRUE.,flag)
  IF(flag > 0) RETURN
  ! ... open file for restart dump
  ufile = ' '
  CALL gfiles(5, ufile, flag)
  IF(flag > 0) RETURN
  ! ... open file and write header with dimensional units to file Plot_timeseries
  CALL gfiles(8,ufile,flag)
  ! ...                   LITHOSTATIC PRESSURE
  ! ... calculate lithostatic pressure as a check on Pressures
  CALL preslith
  ! ... write lithostatic Pressure
  WRITE(fupdef,9095) 'Lithostatic Pressure Limit', unitlabel(11)
  CALL printar(2,plith,ibc,fupdef,unitfac(11),24,000)
  ! ... Flag active seepage nodes based on node pressure
  DO  j = 1,ny
     DO  k = 1,nz
        DO  i = 1,nx
           ik = (k-1)*nx + i
           ! ... Select seepage nodes or combination seepage/precipitation nodes
           IF (ibc(ik,j) == 30 .OR. ibc(ik,j) == 50) THEN
              ic = (k-1)*nx + i + (j-1)*nxz
              IF(p(ic) >= patm) THEN
                 p(ic) = patm
                 seeping(ik,j) = .TRUE.
                 ! ... For now set the temperature to i.c. value
                 ! ... ***  Already done when i.c. data read in
                 !...                    t(ic) =
              ELSE IF(p(ic) < patm) THEN
                 seeping(ik,j) = .FALSE.
              END IF
           END IF
        END DO
     END DO
  END DO
!!$  !..c.....Flag active seepage nodes based on fluid flux
!!$  ! ... Not necessary for i.c.
!!$  !..        DO 500 j = 1,Ny
!!$  !..          DO 490 k = 1,Nz
!!$  !..            DO 480 i = 1,Nx
!!$  !..               ij = (j-1)*Nx + i
!!$  !..               ik = (k-1)*Nx + i
!!$  !..c.....Select  seepage nodes
!!$  !..              IF (ibc(ik,j).EQ.30) THEN
!!$  !..                 ic = (k-1)*Nx + i + (j-1)*Nxz
!!$  !..                 CALL GOVEQN(gveqnm,gveqne,i,j,k)
!!$  !..C---add in the capacitance and source/sink terms
!!$  !..              gveqnm = gveqnm - (Xm(ic)-Xmoldt(ic))/Delt + Q(ic)
!!$  !..              gveqne = gveqne - (En(ic)-Enoldt(ic))/Delt + Qh(ic)
!!$  !..c.....Convert to flux and get correct sign (out is negative)
!!$  !..                qfluxsp = -gveqnm/land_surf_area(ij)
!!$  !..                qhfluxsp = -gveqne/land_surf_area(ij)
!!$  !..                IF (qfluxsp.lt.0._KDP) THEN
!!$  !..                   seeping(ik,j) = .true.
!!$  !..                elseif (qfluxsp.ge.0._KDP) THEN
!!$  !..                   seeping(ik,j) = .false.
!!$  !..                endif
!!$  !..             endif
!!$  !..  480       CONTINUE
!!$  !..  490     CONTINUE
!!$  !..  500   CONTINUE
  ! ... read header record for TIME PERIOD (simulation period)
  IF (informt == 1) READ(fuins,9000)
  RETURN
380 WRITE(fustdout,9165)
    WRITE(fuclog,9165)
9165 FORMAT (/5X, '***** Problem with the variable "newfmt" *****',  &
             /11X, 'Use 0 instead of F, and 1 instead of T')
  ierr(36) = .TRUE.
  RETURN
390 WRITE(fustdout,9170)
  WRITE(fuclog,9170)
9170 FORMAT (/5X, '***** Problem with the DIMENSIONS input line *****'/  &
       11X, 'Check that NO value is assigned to bandwidth')
  ierr(37) = .TRUE.
  RETURN
9000 FORMAT (1X)
!!$9025 FORMAT (/, '--- Time Step Control ---', /,TR10,  &
!!$       'Maximum number of time steps:                      ',i5/TR10,  &
!!$       'Max factor for increasing time step:               ',f5.1/TR10,  &
!!$       'Max % change in nodal pressure per time step:      ',f5.1/TR10,  &
!!$       'Max % change in nodal enthalpy per time step:      ',f5.1/TR10,  &
!!$       'Max change in nodal water satn per time step:      ',f5.3)
!!$  ! ...   , /,TR10,
!!$  ! ... &        'Max no. of cells fluid may traverse per time step   = ',
!!$  ! ... &        G10.3, /, 20X, 'ABOVE OPTION NOT YET IMPLEMENTED')
9040 FORMAT (1P,TR10, 'Max number of SSOR iterations =',i13, /,TR10,  &
       'SSOR tolerance =', 18X, d10.3, /,TR10, 'SSOR relaxation parameter =', d17.3)
9045 FORMAT (10X,  &
       'Potential difference threshold for upstream weighting: ', 1PG9.2/10X,  &
       'Weighting Factor for Temporal Differencing:            ', 0PF9.1)
  !...     &        /, 15X, 'ABOVE OPTION NOT YET IMPLEMENTED, Theta=1 used')
9050 FORMAT (10X, 'Residual Water Saturation:                       ',f8.3/  &
       TR10, 'Residual Steam Saturation:                       ',f8.3)
9055 FORMAT (tr10,a,'Spatially variable (defined below)')
9060 FORMAT (tr10,a,1pg10.3,' - Spatially uniform')
9065 FORMAT (tr10,a,1pg10.3)
9070 FORMAT (/, '--- Radial Coordinates Selected ---', /,  &
       TR10, 1P, 'Well radius = ', g10.3, ' (cm)', /,TR10,  &
       'Maximum reservoir radius = ', g10.3, ' (cm)')
!!$ 9080 FORMAT (1P, 'Row',i4, e13.3, e16.5, e13.4, e16.5, e13.4)
9095 FORMAT (/, '--- ', a, ' ---  ', a)
9130 FORMAT (10X, 'Uniform Value -  ', 1P, g11.4)
9175 FORMAT (10X, a)
9180 FORMAT (10X, 2A)
9190 FORMAT (a)
9195 FORMAT (/, 7X, a)
9245 FORMAT (10X, a, ' Weighting used for Enthalpy & Relative Permeabilty')
9255 FORMAT (10X, a, ' Relative Permeabilty Function Selected')
END SUBROUTINE gdata
