SUBROUTINE gmres(n,msdr,rhs,sol,stop_tol,maxits,aa,ja,ia,alu,jlu,ju,iierr,n_iter,r_norm)
  !                 *** ILU - Preconditioned GMRES ***                  
  !----------------------------------------------------------------------
  ! This is a simple version of the ILU preconditioned GMRES algorithm. 
  ! GMRES uses the L and U matrices generated 
  ! from the subroutine ILU to precondition the GMRES algorithm.        
  ! The preconditioning is applied to the right. The stopping criterion  
  ! utilized is based simply on reducing the relative residual norm to stop_tol
  !     absolute stop_tol (eps_a) is used to handle small initial rhs.
  !                                                                      
  ! USAGE: first call ILUT or ILUK to set up preconditioner and 
  !    then call gmres.                                                    
  !----------------------------------------------------------------------
  ! adapted from  Y. Saad - 5/90
  ! see also chp.9.3.2 of Saad (2003) book for algorithm 9.5
  !----------------------------------------------------------------------
  ! subroutines called :                                           
  ! amux   : SPARSKIT routine to do the matrix*vector multiplication 
  ! lusol : combined forward and backward solves from preconditioning
  ! several BLAS1 routines                                                
  !----------------------------------------------------------------------
  USE machine_constants, ONLY: kdp, epsmac
  USE f_units, ONLY: fuclog, fustdout
  USE control, ONLY: ierr, ioptpr
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n        ! ... The row dimension of A
  INTEGER, INTENT(IN) :: msdr     ! ... Size of Krylov subspace
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: rhs     ! ... right hand side vector
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN OUT) :: sol     ! ... solution vector, initial guess
                                                          ! ...   on input
  REAL(KIND=kdp), INTENT(IN) :: stop_tol    ! ... tolerance for stopping criterion. Iterations stop
  ! ...   when L2-norm(current residual)/L2-norm(initial residual) <= stop_tol
  ! ...   a small absolute tolerance is also used 
  INTEGER, INTENT(IN) :: maxits     ! ... maximum number of iterations allowed
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: aa    ! ... Matrix A in compressed
  INTEGER, DIMENSION(:), INTENT(IN) :: ja           ! ... sparse row format
  INTEGER, DIMENSION(:), INTENT(IN) :: ia
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: alu     ! ...  LU matrix stored in Modified 
                                                      ! ...    Sparse Row (MSR) format
                                                      ! ...  The diagonal (stored in alu(1:n))
                                                      ! ...   is inverted. 
  INTEGER, DIMENSION(:), INTENT(IN) :: jlu    ! ... Each i-th row of alu,jlu matrices
  ! ...   contains the i-th row of L (excluding 
  ! ...   the unit diagonal) followed by the i-th
  ! ...   row of U.
  INTEGER, DIMENSION(:), INTENT(IN) :: ju     ! ... The pointers to the beginning of each row 
  ! ...   of U in the matrices alu,jlu
  INTEGER, INTENT(OUT) :: iierr     ! ... Error message flag
  !            0 : successful convergence to solution
  !            1 : iteration limit reached without convergence
  !           -1 : initial solution gives residual of zero
  INTEGER, INTENT(OUT), OPTIONAL :: n_iter           ! ... Total iterations at convergence
  REAL(kind=kdp), INTENT(OUT), OPTIONAL :: r_norm    ! ... Norm of relative residual at convergence
  !
  ! ... Patch since automatic arrays clash with Java using Visual Fortran90 v6.0
!!$  REAL(KIND=kdp), DIMENSION(n,msdr+1), TARGET :: vv     ! ... work array. stores the Arnoldi
  REAL(KIND=kdp), DIMENSION(:,:), ALLOCATABLE, TARGET :: vv    ! ... work array. stores the Arnoldi
                                                               ! ...   basis vectors
  !
  REAL(KIND=kdp), DIMENSION(:), POINTER :: vvp1   ! ... work array pointer to vector slice
  REAL(KIND=kdp), DIMENSION(:), POINTER :: vvp2   ! ... work array pointer to vector slice
  REAL(KIND=kdp), EXTERNAL :: ddot, dnrm2
  !
  INTEGER :: i, i1, ii, itno, j, jj, k, k1
  INTEGER :: a_err, da_err
  REAL(KIND=kdp) :: eps1, gam, ro, t
!!$  REAL(KIND=kdp), DIMENSION(msdr+1,msdr) :: hh
!!$  REAL(KIND=kdp), DIMENSION(msdr) :: c, s
!!$  REAL(KIND=kdp), DIMENSION(msdr+1) :: rs
  REAL(KIND=kdp), DIMENSION(:,:), ALLOCATABLE :: hh
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: c, s, rs, mvy
  REAL(KIND=kdp), PARAMETER :: eps_a=1.e-16_kdp   ! *** make multiple of mach eps
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$RCSfile: gmres.f90,v $//$Revision: 1.12 $'
  !     ------------------------------------------------------------------
  !...
  ! ... comments follow Templates p.20 as closely as possible
  ! ... see also algorithm 9.5 of Saad (2003) book
  ! ... allocate work space
  ALLOCATE(vv(n,msdr+1), hh(msdr+1,msdr), c(msdr), s(msdr), rs(msdr+1), mvy(n),  &
       STAT=a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: gmres"
     ierr(196) = .TRUE.
     RETURN
  ENDIF
  ! ...  compute initial residual vector; inside the outer loop in Templates
  vvp1 => vv(:,1)
  CALL amux(sol, vvp1, aa, ja, ia)
  DO  j=1,n
     vv(j,1) = rhs(j) - vv(j,1)      ! ... initial residual vector
  END DO
  itno = 0                           ! ... iteration counter
  DO           ! ... outer loop of iteration
     vvp1 => vv(:,1)
     ro = dnrm2(n, vvp1, 1)         ! ... current residual norm
     IF(itno == 0) THEN
        IF(ioptpr(3) == 2) THEN
           WRITE(fuclog,2001) 'GMRES:','Iteration ',itno,  &
                '; Norm of residual .... ',ro
           WRITE(fustdout,2001) 'GMRES:','Iteration ',itno,  &
                '; Norm of residual .... ',ro
2001       FORMAT(tr5,a/tr5,a,i4,a,1pe12.5)
        END IF
        IF(PRESENT(n_iter)) n_iter = itno
        IF(PRESENT(r_norm)) r_norm = 1._kdp
        eps1 = stop_tol*ro     ! ... for later convergence test
     END IF
     ! ... Drop this test. At steady state, null solution vector is expected.
!!$     IF(ro <= stop_tol + eps_a) THEN     ! ... residual less than tolerance, sol is sufficient
!!$        iierr = -1
!!$        RETURN
!!$     END IF
     t = 1.0_kdp/ro
     DO  j=1,n
        vv(j,1) = vv(j,1)*t       ! ...  v-1 vector; initial vector of Krylov space
     END DO
     ! ... initialize first term of rhs of hessenberg system; s in Templates
     rs(1) = ro
     DO i=1,msdr          ! ... inner loop of iteration; msdr steps between restarts
        itno = itno + 1
        i1 = i + 1
        vvp1 => vv(:,i)
        vvp2 => vv(:,i1)
        CALL lusol(vvp1, rhs, alu, jlu, ju)
        CALL amux(rhs, vvp2, aa, ja, ia)         ! ... vvp2 contains w of Templates
        ! ... Build the orthogonal Krylov vector space using modified Gram-Schmidt algorithm
        ! ...     applying the Arnoldi algorithm
        DO  j=1,i
           vvp1 => vv(:,j)
           t = ddot(n, vvp1, 1, vvp2, 1)
           hh(j,i) = t
           CALL daxpy(n, -t, vvp1, 1, vvp2, 1)
        END DO
        t = dnrm2(n, vvp2, 1)
        hh(i1,i) = t
        IF (t /= 0.0_kdp) THEN        ! ... test for breakdown but no reothorgonalization done
           t = 1.0_kdp/t
           DO   k=1,n
              vv(k,i1) = vv(k,i1)*t
           END DO
        END IF
        NULLIFY(vvp1,vvp2)
        ! ... vv vectors contain the set of orthonormal basis vectors
        !     done with modified gram schimdt and arnoldi step
        !     now  update QR factorization of hh, Hessenberg matrix,
        ! ...    using Givens rotations
        ! ... Apply previous rotations to i-th column of hh
        DO  k=2,i
           k1 = k-1
           t = hh(k1,i)
           hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
           hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
        END DO
        gam = SQRT(hh(i,i)**2 + hh(i1,i)**2)
        !     if gamma is zero then any small value will do
        !     will affect only residual estimate; seems like only hh(i,i) may be zero
        ! ...       by Greenbaum p.40. maybe breakdown makes gamma be zero.
        IF (gam <= 0.0_kdp) gam = epsmac
        ! ... calculate next rotation
        c(i) = hh(i,i)/gam       ! ... Kelley algorithm
        s(i) = hh(i1,i)/gam      ! ... more robust than Greenbaum
        rs(i1) = -s(i)*rs(i)     ! ... s in Templates
        rs(i) =  c(i)*rs(i)
        hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)   ! ... in Greenbaum and Kelley not Templates
        ! ... Calculate residual norm and test for convergence
        ro = ABS(rs(i1))
        IF(ioptpr(3) == 2) THEN
           WRITE(fuclog,2002) 'Iteration ',itno,  &
                '; Norm of residual .... ',ro
           WRITE(fustdout,2002) 'Iteration ',itno,  &
                '; Norm of residual .... ',ro
2002       FORMAT(tr5,a,i4,a,1pe12.5)
        END IF
        IF(PRESENT(n_iter)) n_iter = itno
        IF(PRESENT(r_norm)) r_norm = (ro*stop_tol)/eps1
        IF (ro <= eps1 + eps_a) EXIT     ! ... convergence on relative residual 
     END DO
     i = MIN(i,msdr)
     ! ... Update step in Templates, for either i or msdr step
     ! ... First solve upper triangular Hessenberg system
     rs(i) = rs(i)/hh(i,i)
     DO  ii=2,i
        k = i-ii+1
        k1 = k+1
        t=rs(k)
        DO  j=k1,i
           t = t - hh(k,j)*rs(j)
        END DO
        rs(k) = t/hh(k,k)     ! ... rs now contains y of Templates
     END DO
     !     form linear combination of vv(*,i)'s to get updated solution
     t = rs(1)
     DO  k=1,n
        rhs(k) = t*vv(k,1)
     END DO
     DO  j=2,i
        t = rs(j)
        DO  k=1,n
           rhs(k) = rhs(k) + t*vv(k,j)     ! ... rhs contains solution to precond system
        END DO                             ! ...   without the initial guess
     END DO
     ! ... apply preconditioning operation to recover solution to original 
     ! ...     unpreconditioned problem without initial guess
     CALL lusol(rhs, mvy, alu, jlu, ju)
     ! ... add in initial guess
     DO  k=1,n
        sol(k) = sol(k) + mvy(k)     ! ... current solution
     END DO
     rhs = mvy
     IF (ro <= eps1 + eps_a) EXIT          ! ... convergence on relative residual
     IF (itno >= maxits) THEN              ! ... itno is incremented in multiples of msdr
        iierr = 1     ! ... iteration limit reached
        RETURN
     END IF
     ! ... Compute residual vector and continue; done before convergence test in Templates
     ! ... Uses results of minimization not explicit residual formula
     ! ... Seems redundant
     DO  j=1,i
        jj = i1-j+1
        rs(jj-1) = -s(jj-1)*rs(jj)
        rs(jj) = c(jj-1)*rs(jj)
     END DO
     vvp2 => vv(:,1)
     DO  j=1,i1
        t = rs(j)
        IF (j == 1)  t = t-1.0_kdp
        vvp1 => vv(:,j)
        CALL daxpy(n, t, vvp1, 1, vvp2, 1)
     END DO
  END DO                           ! ... end outer loop
  iierr = 0     ! ... convergence
  DEALLOCATE(vv, hh, c, s, rs, mvy,  &
       STAT=da_err)
  IF (da_err /= 0) THEN  
     PRINT *, "Array deallocation failed: gmres"
     ierr(191) = .TRUE.
  ENDIF

CONTAINS

  SUBROUTINE lusol(y, x, alu, jlu, ju)
    ! ... Solves the system LU*x = y,
    ! given an LU decomposition of a matrix stored in (alu, jlu, ju)
    ! in modified sparse row format MSR
    IMPLICIT NONE
    REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: y     ! ... right hand side vector
    REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: x    ! ... solution vector
    REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: alu   ! ... LU matrix stored in Modified 
                                                      ! ...   Sparse Row (MSR) format
                                                      ! ... Provided by ILU routine
    INTEGER, DIMENSION(:), INTENT(IN) :: jlu    ! ... Each i-th row of alu,jlu matrices
    ! ...   contains the i-th row of L (excluding 
    ! ...   the unit diagonal) followed by the i-th
    ! ...   row of U.
    INTEGER, DIMENSION(:), INTENT(IN) :: ju     ! ... The pointers to the beginning of each row 
                                                ! ...   of U in the matrices alu,jlu
    !
    INTEGER :: i, k
    ! ... Set string for use with RCS ident command
    !CHARACTER(LEN=80) :: ident_string='$RCSfile: gmres.f90,v $//$Revision: 1.12 $'
    !     ------------------------------------------------------------------
    !...
    ! ... forward solve
    DO  i=1,n              ! ... n is known from host
       x(i) = y(i)
       DO  k=jlu(i),ju(i)-1
          x(i) = x(i) - alu(k)*x(jlu(k))
       END DO
    END DO
    ! ... backward solve
    DO  i=n,1,-1
       DO  k=ju(i),jlu(i+1)-1
          x(i) = x(i) - alu(k)*x(jlu(k))
       END DO
       x(i) = alu(i)*x(i)
    END DO
  END SUBROUTINE lusol

  SUBROUTINE amux(x,y,a,ja,ia)
    ! ... Calculates product of matrix A times vector x
    ! multiplies a matrix by a vector using the dot product form
    ! Matrix A is stored in compressed sparse row storage, CSR.
    ! ... y = A*x
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: x     ! ... vector x
    REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: y    ! ... result vector y 
    REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: a     ! ... Matrix A in compressed
    INTEGER, DIMENSION(:), INTENT(IN) :: ja           ! ... sparse row format
    INTEGER, DIMENSION(:), INTENT(IN) :: ia
    !
    INTEGER :: i, k
    REAL(KIND=kdp) :: t
    ! ... Set string for use with RCS ident command
    !CHARACTER(LEN=80) :: ident_string='$RCSfile: gmres.f90,v $//$Revision: 1.12 $'
    !     ------------------------------------------------------------------
    !...
    DO  i = 1,n                ! ... n is known from host
       ! ... Compute the inner product of row i with vector x
       t = 0.0_kdp
       DO  k=ia(i),ia(i+1)-1
          t = t + a(k)*x(ja(k))
       END DO
       y(i) = t
    END DO
  END SUBROUTINE amux

END SUBROUTINE gmres
