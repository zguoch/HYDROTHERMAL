SUBROUTINE gmres_asem_solve
  ! ... Purpose:  Assemble the finite-difference equation matrices for the mesh
  ! ...           and then solve them using a preconditioned Generalized Minimum
  ! ...           Residual Method
  ! ...           The solution is returned in xxs
  USE machine_constants, ONLY: kdp
  USE f_units, ONLY: fuclog, fustdout
  USE control, ONLY: ierr, ioptpr, itmcntl, tmcut, tot_iter_gm, iter_gm_max
  USE ilupc_mod
  USE mesh
  USE solver_gmres
  USE variables
  IMPLICIT NONE
  INTERFACE
     SUBROUTINE rowscale(nrow,job,norm,a,ja,ia,diag,b,jb,ib,ierr)
       USE machine_constants, ONLY: kdp
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: nrow
       INTEGER, INTENT(IN) :: job
       INTEGER, INTENT(IN) :: norm
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: a
       INTEGER, DIMENSION(:), INTENT(IN) :: ja
       INTEGER, DIMENSION(:), INTENT(IN) :: ia
       REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: diag
       REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: b
       INTEGER, DIMENSION(:), INTENT(OUT) :: jb
       INTEGER, DIMENSION(:), INTENT(OUT) :: ib
       INTEGER, INTENT(OUT) :: ierr
     END SUBROUTINE rowscale
  END INTERFACE     !!$*** for flint
  INTERFACE         !!$*** for flint
     SUBROUTINE colscale(nrow,job,norm,a,ja,ia,diag,b,jb,ib,ierr)
       USE machine_constants, ONLY: kdp
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: nrow    ! ... The row dimension of A
       INTEGER, INTENT(IN) :: job     ! ... job selector. 0: calculate B only
       INTEGER, INTENT(IN) :: norm     ! ... norm selector. 1: 1-norm,
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: a     ! ... Matrix A in compressed
       INTEGER, DIMENSION(:), INTENT(IN) :: ja           ! ... sparse row format
       INTEGER, DIMENSION(:), INTENT(IN) :: ia
       REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: diag  ! ... Diagonal matrix of scale factors
       REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: b        ! ... Matrix B in compressed
       INTEGER, DIMENSION(:), INTENT(OUT) :: jb              ! ... sparse row format
       INTEGER, DIMENSION(:), INTENT(OUT) :: ib
       INTEGER, INTENT(OUT) :: ierr              ! ... error message. 0 : Normal return
     END SUBROUTINE colscale
  END INTERFACE     !!$*** for flint
  INTERFACE         !!$*** for flint
     SUBROUTINE gmres(n,msdr,rhs,sol,stop_tol,maxits,aa,ja,ia,alu,jlu,ju,ierr,  &
          n_iter,r_norm)
       USE machine_constants, ONLY: kdp
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n
       INTEGER, INTENT(IN) :: msdr
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: rhs
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: sol
       REAL(KIND=kdp), INTENT(IN) :: stop_tol
       INTEGER, INTENT(IN) :: maxits
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: aa
       INTEGER, DIMENSION(:), INTENT(IN) :: ja
       INTEGER, DIMENSION(:), INTENT(IN) :: ia
       REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: alu
       INTEGER, DIMENSION(:), INTENT(IN) :: jlu
       INTEGER, DIMENSION(:), INTENT(IN) :: ju
       INTEGER, INTENT(OUT) :: ierr
       INTEGER, INTENT(OUT), OPTIONAL :: n_iter
       REAL(KIND=kdp), INTENT(OUT), OPTIONAL :: r_norm
     END SUBROUTINE gmres
  END INTERFACE
  !
  INTEGER :: i, ic, iierr, job, ma, n_iter, norm
  INTEGER :: a_err, da_err
  REAL(KIND=kdp) :: corvech, corvecp
  ! ... Patch since automatic arrays clash with Java using Visual Fortran90 v6.0
!!$  REAL(KIND=kdp), DIMENSION(nxyz) :: dho, dpo
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: dho, dpo
!!$  REAL(KIND=kdp), DIMENSION(2*mrnomax) :: urhs
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: urhs
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.5 $//$Date: 2010/03/30 23:57:14 $'
  !     ------------------------------------------------------------------
  ! ... allocate local storage
  ALLOCATE (dho(nxyz), dpo(nxyz), urhs(2*mrnomax),  &
       STAT=a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: gmres_asem_solve"
     ierr(195) = .TRUE.
     RETURN
  ENDIF
  ! ... initialize arrays for change in pressure and enthalpy (new and old) 
  dpo = 0._kdp
  dpn = 0._kdp
  dho = 0._kdp
  dhn = 0._kdp
  ! ... Initialize the variables to record the residual P and H from N-R solution??
  corvecp = 0._kdp
  corvech = 0._kdp
  ! ... Assemble the N-R difference equations
  CALL assemble_nre_csr
  ! ... Scale the matrix equations
  job = 1
  norm = 2
  neq = 2*mrnomax
  CALL rowscale(neq,job,norm,av,ja,ia,diagr,bv,jb,ib,iierr)
  IF(iierr /= 0) THEN
     WRITE(fuclog,*) '*** ERROR: Scaling row no. ', iierr
     ierr(81) = .TRUE.
!!$     RETURN
     go to 999
  END IF
  CALL colscale(neq,job,norm,bv,jb,ib,diagc,av,ja,ia,iierr)
  IF(iierr /= 0) THEN
     WRITE(fuclog,*) '*** ERROR: Scaling column no. ', iierr
     ierr(81) = .TRUE.
!!$     RETURN
     go to 999
  END IF
  DO i=1,neq
     urhs(i) = diagr(i)*rhs(i)
  END DO
  ! ... Solve linear N-R matrix equation
  xxs = 0._kdp     ! ... Initial guess for solution vector
  SELECT CASE (ilu_method)
  CASE (1)
     CALL ilut(neq,av,ja,ia,lev_fill,drop_tol,ailu,jailu,iilu,iierr)
     IF(ioptpr(3) == 2) THEN
        WRITE(fuclog,9001) 'Number of non-zero elements in A .......... ',  &
             ia(neq+1) - ia(1)
        WRITE(fuclog,9001) 'Number of non-zero elements in ILU(tol) ... ',  &
             jailu(neq+1) - jailu(1) + neq
9001    FORMAT(tr5,a,i8)
     END IF
  CASE (2)
     CALL iluk(neq,av,ja,ia,lev_fill,ailu,jailu,iilu,iierr)
  END SELECT
  IF(iierr /= 0) THEN
     WRITE(fuclog,*) 'Error in Preconditioning: ', iierr
     SELECT CASE (iierr)
     CASE (1:)
        ierr(87) = .TRUE.
     CASE (-1)
        ierr(82) = .TRUE.
     CASE (-2)
        ierr(83) = .TRUE.
     CASE (-3)
        ierr(84) = .TRUE.
     CASE (-4)
        ierr(85) = .TRUE.
     CASE (-5)
        ierr(86) = .TRUE.
     CASE (-6)
        ierr(194) = .TRUE.
     CASE (-7)
        ierr(193) = .TRUE.
     END SELECT
!!$     RETURN
     go to 999
  END IF
  !  urhs = rhs
  CALL gmres(neq,m_save_dir,urhs,xxs,stop_tol_gmres,maxit_gmres,av,ja,ia,  &
       ailu,jailu,iilu,iierr,  &
       n_iter,r_norm_gmres)
  WRITE(fustdout,9005) 'GMRES solver return status: ', iierr,  &
       '; Number of iterations ... ', n_iter,  &
       'Relative norm of residual ............................. ', r_norm_gmres
  WRITE(fuclog,9005) 'GMRES solver return status: ', iierr,  &
       '; Number of iterations ... ', n_iter,  &
       'Relative norm of residual ............................. ', r_norm_gmres
9005 FORMAT(tr5,a,i1,a,i4/tr5,a,1pg11.4)
  tot_iter_gm = tot_iter_gm + n_iter
  iter_gm_max = MAX(iter_gm_max,n_iter)
  IF(iierr /= 0) THEN
     WRITE(fuclog,*) 'Error in GMRES solving: ', iierr
!     SELECT CASE (iierr)
!     CASE (1)                !  This should not be counted as a run-stopping error. N-R iteration will restart with a smaller time step size.
!        ierr(88) = .TRUE.
!     CASE (-1)               !  This case has been disabled in gmres
!        ierr(89) = .TRUE.
!     END SELECT
!!$     RETURN
!   go to 999   ! continue to cut time step
  END IF
  ! ... Cut time step if max number of GMRES iterations allowed is reached
  ! ...       without convergence
  IF (n_iter >= maxit_gmres) THEN
     IF (ioptpr(3) >= 1) THEN
        WRITE (fustdout,'(2A)') '*** MAXIMUM NUMBER of GMRES ITERATIONS ',  & 
             'REACHED WITHOUT CONVERGENCE, time step cut ***'
        WRITE (fuclog,'(2A)') '*** MAXIMUM NUMBER of GMRES ITERATIONS ',  &
             'REACHED WITHOUT CONVERGENCE, time step cut ***'
     END IF
     tmcut = .TRUE.
     itmcntl(21) = itmcntl(21) + 1
!!$     RETURN
     go to 999
  END IF
  ! ... Extract the solutions from the xxs vector and descale
  DO ic=1,nxyz
     ma = mrno(ic)
     IF(ma > 0) THEN
        dpn(ic) = diagc(2*ma-1)*xxs(2*ma-1)
        dhn(ic) = diagc(2*ma)*xxs(2*ma)
     END IF
  END DO
999 continue
  DEALLOCATE(dho, dpo, urhs,  &
       STAT=da_err)
  IF (da_err /= 0) THEN  
     PRINT *, "Array deallocation failed: gmres_asem_solve"
     ierr(192) = .TRUE.
  ENDIF
END SUBROUTINE gmres_asem_solve
