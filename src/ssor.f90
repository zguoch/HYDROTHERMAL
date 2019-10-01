SUBROUTINE ssor
  ! ... Purpose:  To assemble N-R matrix equations for each
  ! ...           X-Z slice and then to solve them using Slice
  ! ...           Successive Over-Relaxation
  ! ...             SSOR first calls FORMEQ which assembles the coefficient
  ! ...           matrix and the right hand side that doesn't change
  ! ...           for each slice J.  It then calls SOLVE to triangularize
  ! ...           the matrix for each slice.  The column vector R is modified
  ! ...           to include Y term residuals and the partials of the Mass
  ! ...           and Energy govn eqns w.r.t. P and H.  Each matrix is then
  ! ...           solved by back substitution by first assuming no flow in
  ! ...           the Y direction.  The solution is returned in R
  USE machine_constants, ONLY: kdp
  USE f_units
  USE parameters          !!$*** for flint
  USE control
  USE fdeq
  USE mesh
!!$  USE parameters
  USE solver
  USE variables
  IMPLICIT NONE
  INTEGER :: i, ic, icym, icyp, j, k, lssor, ml, mm, mym, myp, nbby
  INTEGER :: a_err, da_err, ii, nn
  REAL(KIND=kdp) :: corvech, corvecp, difmx, dipmx, wosave
  ! ... Patch since automatic arrays clash with Java using Visual Fortran90 v6.0
!!$  REAL(KIND=kdp), DIMENSION(2*nxz,mbw,ny) :: abig
!!$  REAL(KIND=kdp), DIMENSION(2*nxz,ny) :: rbig
!!$  REAL(KIND=kdp), DIMENSION(nxyz) :: dho, dpo
  REAL(KIND=kdp), DIMENSION(:,:,:), ALLOCATABLE :: abig
  REAL(KIND=kdp), DIMENSION(:,:), ALLOCATABLE :: rbig
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: dho, dpo
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.8 $//$Date: 2007/10/19 17:54:31 $'
  !     ------------------------------------------------------------------
  ! ... allocate local storage
  ALLOCATE (abig(2*nxz,mbw,ny), rbig(2*nxz,ny), dho(nxyz), dpo(nxyz),  &
             STAT=a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: ssor"
     ierr(197) = .TRUE.
     RETURN
  ENDIF
  ! ... initialize change in pressure and enthalpy (new and old)
  ! ... this assumes that there is no flow in the Y direction on the first iteration
  dpo = 0._kdp
  dpn = 0._kdp
  dho = 0._kdp
  dhn = 0._kdp
  IF (ny == 1) THEN       ! ...  Only one slice,  skip SSOR
     ! ... ---------- Matrix Solution for Only One Slice ---------------
     ! ...   initialize variables to hold the delta P and delta H from N-R solution
     corvecp = 0._kdp
     corvech = 0._kdp
     j = 1
     CALL formeq(j)
     nbby = 2*nb(j)
     ! ... Factor matrix A and solve using direct band matrix solver
     CALL solve(3,nbby)
     ! ... Descale the solution
     nn = nbby
     DO  ii = 1,nn
        rlit(ii) = cscal(ii)*rlit(ii)
     END DO
     ! ... Load the solution back into the dependent variables
     DO  ml = 1,nb(j)
        i = npp(2*ml-1,j)
        k = npp(2*ml,j)
        ic = (k-1)*nx + i + (j-1)*nxz
        dpn(ic) = rlit(2*ml-1)
        dhn(ic) = rlit(2*ml)
        ! ...   IF (ABS(CORVECP).LT.ABS(RLIT(2*ML-1))) CORVECP=RLIT(2*ML-1)
        ! ...   IF (ABS(CORVECH).LT.ABS(RLIT(2*ML))) CORVECH=RLIT(2*ML)
     END DO
!!$!***special output
!!$      open(61,file='x.out')
!!$      write(61,'(1pg25.15)')  (rlit(i),i=1,2*nb(j))
!!$      close(61)
!!$!***end special output
  ELSEIF (ny > 1) THEN
! ***** this ssor is broken until the proper row and column scaling is done
! *****      and both equations are used in the convergence test
     ! ... ---------------- Beginning of SSOR Iteration ---------------------
     ! ...    Set up a matrix for each X-Z slice and triangularize it
     DO  j = 1,ny
        ! ... Skip slice if there are no active nodes (ie. all constant value nodes)
        IF (nb(j) == 0) CYCLE
        ! ...      set up coefficient matrix and right hand side matrix
        ! ...         that doesn't change for this slice, J
        CALL formeq(j)
        ! ...    store right hand side matrix for this slice, J, in the Big RHS
        ! ...        array (RBIG); the number of elements in the matrix R is NBBY
        ! ...        which is equal to 2x the # of active nodes in slice, J: NB(J)
        nbby = 2*nb(j)
        DO  ml = 1,nbby
           rbig(ml,j) = rlit(ml)
        END DO
        ! ...    triangularize matrix A for slice J
        CALL solve(1,nbby)
        ! ...    store the matrix for slice J in the big A array (ABIG)
        nbby = 2*nb(j)
        DO  ml = 1,nbby
           DO  mm = 1,mbw
              abig(ml,mm,j) = alit(ml,mm)
           END DO
        END DO
     END DO
     ! ... -------------------- End of matrix set up ----------------------
     ! ...         Begin SSOR iteration - iterate until the max number of
     ! ...               iterations per timestep, or until the maximum
     ! ...               residual is less than the tolerance
     ! ...         overrelaxation factor - WO & WOSAVE
     wosave = wo
     DO  lssor = 1,lssormax
        IF (lssor == 1) THEN
           wo = 1._kdp
        ELSE
           wo = wosave
        END IF
        ! ...       initialize variables to record maximum change of P into Y???
        dipmx = 0._kdp
        difmx = 0._kdp
        ! ...   initialize variables to hold the delta P and delta H from N-R solution
        corvecp = 0._kdp
        corvech = 0._kdp
        ! ...     Modify RHS column matrix for each slice J by adding in the
        ! ...     Y direction terms and solve by back substitution for each slice
        DO  j = 1,ny
           ! ...      skip slice if there are no active nodes (ie. all constant value nodes)
           IF (nb(j) == 0) CYCLE
           DO  ml = 1,nb(j)
              i = npp(2*ml-1,j)
              k = npp(2*ml,j)
              ic = (k-1)*nx + i + (j-1)*nx*nz
              icyp = ic + nxz
              icym = ic - nxz
              mym = ic
              myp = mym + nxz
              ! ...    pass back the RHS column matrix R for this slice from RBIG
              rlit(2*ml) = rbig(2*ml,j)
              rlit(2*ml-1) = rbig(2*ml-1,j)
              IF (j /= ny) THEN
                 ! ...         governing equations
                 rlit(2*ml-1) = rlit(2*ml-1) - ty(myp)*  &
                      ((rs(icyp)*usy(myp)+rs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                      (rw(icyp)*uwy(myp)+rw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*  &
                      (p(icyp)-p(ic))
                 rlit(2*ml) = rlit(2*ml) - ty(myp)*  &
                      ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                      (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*  &
                      (p(icyp)-p(ic)) - tyk(myp)*(tc(icyp)-tc(ic))
                 ! ...      derivatives of mass eqn w.r.t. P & H
                 rlit(2*ml-1) = rlit(2*ml-1) - dpn(icyp)*ty(myp)*  &
                      ((rs(icyp)*usy(myp)+rs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                      (rw(icyp)*uwy(myp)+rw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                      ((dyvdsp(myp)*rs(icyp)+drsp(icyp)*yvds(myp))*usy(myp)+  &
                      dyvdsp(myp)*rs(ic)*(1._kdp-usy(myp))+  &
                      (dyvdwp(myp)*rw(icyp)+drwp(icyp)*yvdw(myp))*uwy(myp)+  &
                      dyvdwp(myp)*rw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic)))  &
                      - dhn(icyp)*ty(myp)*((dyvdsh(myp)*rs(icyp)+drsh(icyp)*yvds(myp))*  &
                      usy(myp)+ dyvdsh(myp)*rs(ic)*(1._kdp-usy(myp))+  &
                      (dyvdwh(myp)*rw(icyp)+drwh(icyp)*yvdw(myp))* uwy(myp)+  &
                      dyvdwh(myp)*rw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic))
                 ! ...      derivatives of energy eqn w.r.t. P & H
                 rlit(2*ml) = rlit(2*ml) - dpn(icyp)*ty(myp)*  &
                      ((hrs(icyp)*usy(myp)+hrs(ic) *(1._kdp-usy(myp)))*yvds(myp)+  &
                      (hrw(icyp)*uwy(myp)+hrw(ic) *(1._kdp-uwy(myp)))*yvdw(myp)+  &
                      ((dyvdsp(myp)*hrs(icyp)+dhrsp(icyp)*yvds(myp))* usy(myp)+  &
                      dyvdsp(myp)*hrs(ic)*(1._kdp-usy(myp))+  &
                      (dyvdwp(myp)*hrw(icyp)+dhrwp(icyp)*yvdw(myp))* uwy(myp)+  &
                      dyvdwp(myp)*hrw(ic)*(1._kdp-uwy(myp)))* (p(icyp)-p(ic)))  &
                      - dpn(icyp)*tyk(myp)*dtp(icyp)
                 rlit(2*ml) = rlit(2*ml) - dhn(icyp)*ty(myp)*  &
                      ((dyvdsh(myp)*hrs(icyp)+dhrsh(icyp)*yvds(myp))* usy(myp)+  &
                      dyvdsh(myp)*hrs(ic)*(1._kdp-usy(myp))+  &
                      (dyvdwh(myp)*hrw(icyp)+dhrwh(icyp)*yvdw(myp))* uwy(myp)+  &
                      dyvdwh(myp)*hrw(ic)*(1._kdp-uwy(myp)))* (p(icyp)-p(ic))  &
                      - dhn(icyp)*tyk(myp)* dth(icyp)
              END IF
              IF (j /= 1) THEN
                 ! ...         governing equations
                 rlit(2*ml-1) = rlit(2*ml-1) - ty(mym)*  &
                      ((rs(ic)*usy(mym)+rs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                      (rw(ic)*uwy(mym)+rw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*  &
                      (p(icym)-p(ic))
                 rlit(2*ml) = rlit(2*ml) - ty(mym)*  &
                      ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                      (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*  &
                      (p(icym)-p(ic)) - tyk(mym)*(tc(icym)-tc(ic))
                 ! ...      derivatives of mass eqn w.r.t. P & H
                 rlit(2*ml-1) = rlit(2*ml-1) - dpn(icym)*ty(mym)*  &
                      ((rs(ic)*usy(mym)+rs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                      (rw(ic)*uwy(mym)+rw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                      (dyvdsp(mym)*rs(ic)*usy(mym)+  &
                      (dyvdsp(mym)*rs(icym)+drsp(icym)*yvds(mym))* (1._kdp-usy(mym))+  &
                      dyvdwp(mym)*rw(ic)*uwy(mym)+  &
                      (dyvdwp(mym)*rw(icym)+drwp(icym)*yvdw(mym))* (1._kdp-uwy(mym)))*  &
                      (p(icym)-p(ic))) - dhn(icym)*ty(mym)*  &
                      (dyvdsh(mym)*rs(ic)*usy(mym)+  &
                      (dyvdsh(mym)*rs(icym)+drsh(icym)*yvds(mym))* (1._kdp-usy(mym))+  &
                      dyvdwh(mym)*rw(ic)*uwy(mym)+  &
                      (dyvdwh(mym)*rw(icym)+drwh(icym)*yvdw(mym))* (1._kdp-uwy(mym)))*  &
                      (p(icym)-p(ic))
                 ! ...      derivatives of energy eqn w.r.t. P & H
                 rlit(2*ml) = rlit(2*ml) - dpn(icym)*ty(mym)*  &
                      ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                      (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                      (dyvdsp(mym)*hrs(ic)*usy(mym)+  &
                      (dyvdsp(mym)*hrs(icym)+dhrsp(icym)*yvds(mym))* (1._kdp-usy(mym))+  &
                      dyvdwp(mym)*hrw(ic)*uwy(mym)+  &
                      (dyvdwp(mym)*hrw(icym)+dhrwp(icym)*yvdw(mym))* (1._kdp-uwy(mym)))*  &
                      (p(icym)-p(ic))) - dpn(icym)*tyk(mym)*  &
                      dtp(icym)
                 rlit(2*ml) = rlit(2*ml) - dhn(icym)*ty(mym)*  &
                      (dyvdsh(mym)*hrs(ic)*usy(mym)+  &
                      (dyvdsh(mym)*hrs(icym)+dhrsh(icym)*yvds(mym))* (1._kdp-usy(mym))+  &
                      dyvdwh(mym)*hrw(ic)*uwy(mym)+  &
                      (dyvdwh(mym)*hrw(icym)+dhrwh(icym)*yvdw(mym))* (1._kdp-uwy(mym)))*  &
                      (p(icym)-p(ic)) - dhn(icym)*tyk(mym)*  &
                      dth(icym)
              END IF
              ! ...    pass back the matrix A for this slice from ABIG
              DO  mm = 1,mbw
                 alit(2*ml-1,mm) = abig(2*ml-1,mm,j)
                 alit(ml*2,mm) = abig(2*ml,mm,j)
              END DO
           END DO
           nbby = 2*nb(j)
           ! ...       solve matrix for slice J by updating RHS & using back substitution
           CALL solve(2,nbby)
           DO  ml = 1,nb(j)
              i = npp(2*ml-1,j)
              k = npp(2*ml,j)
              ic = (k-1)*nx + i + (j-1)*nx*nz
              dpo(ic) = dpn(ic)
              dho(ic) = dhn(ic)
              ! ...               residual times relaxation factor
              dpn(ic) = rlit(2*ml-1)*wo + (1._kdp-wo)*dpo(ic)
              dhn(ic) = rlit(2*ml)*wo + (1._kdp-wo)*dho(ic)
              difmx = MAX(difmx,ABS(dpn(ic)-dpo(ic)))
              dipmx = MAX(dipmx,ABS(dpn(ic)))
              corvecp = MAX(corvecp,ABS(rlit(2*ml-1)))
              corvech = MAX(corvech,ABS(rlit(2*ml)))
           END DO
        END DO
        IF (ioptpr(3) >= 1) THEN
           WRITE (fuclog, 9010) ltime, lnri, lssor, difmx/dipmx
9010       FORMAT (1P, i4, i3, i3,  &
                ' - SSOR Iteration,  relative change in P correction = ', g12.5)
           IF (ioptpr(3) == 2) WRITE (fuclog, 9005) ltime, lnri, lssor,  &
                corvecp, corvech
9005       FORMAT (1P, i4, i3, i3, ' Max correction vector:  P = ', g12.5,  &
                ';  H = ', g12.5)
        END IF
        ! ... Convergence test
        IF (difmx/dipmx < tol) exit
     END DO      ! ... End of SSOR Iteration
     if(lssor > lssormax)  &          ! ... Convergence failure
          WRITE (fuclog, 9015) 'Maximum Number,',lssormax,  &
          ', of SSOR Iterations Reached ', 'Without Meeting Tolerance of ', tol
9015 FORMAT (tr2,a,i5,2A,1PG11.4)
  ENDIF
  DEALLOCATE(abig, rbig, dho, dpo,  &
       STAT=da_err)
  IF (da_err /= 0) THEN  
     PRINT *, "Array deallocation failed: ssor"
     ierr(189) = .TRUE.
  ENDIF
END SUBROUTINE ssor
