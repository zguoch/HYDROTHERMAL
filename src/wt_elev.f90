SUBROUTINE wt_elev
  ! ... Computes the elevation of the water table using linear interpolation
  ! ...      of the pressure field to the elevation of atmospheric pressure
  ! ... This is approximate because with cell centered mesh the slope of the 
  ! ...     pressure is discontinuous across the cell boundary
  USE machine_constants, ONLY: kdp
!  USE bc, ONLY: ibc, seeping, cdtn, denflux, ehflux, qprecip, nhcond, nprecip, nseep
!  USE control
  USE mesh
  USE parameters
  USE variables, ONLY: denw, fs_elev, ktop, p, satnw, zz, dfs_l2norm
  USE control, ONLY: ierr
  USE units
  IMPLICIT NONE
  INTEGER :: i, ic, ij, j, k, k0, k1
  REAL(KIND=kdp) :: p0, p1, z0, z1
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: fs_elev1, fs_elev2, del_fs
  INTEGER :: alloc_stat
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5 $//$Date: 2003/10/15 22:01:44 $'
  ALLOCATE(fs_elev1(nxy), fs_elev2(nxy), del_fs(nxy), STAT=alloc_stat)
  If (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: wt_elev'
    ierr(175) = .TRUE.
    RETURN
  END IF
  !     ------------------------------------------------------------------
  !...
  ! ... save the previous water-table (free-surface) elevation
  fs_elev1 = 0._kdp
  fs_elev2 = 0._kdp
  DO  j = 1,ny
     DO  i = 1,nx
        ij = (j-1)*nx + i
        DO  k = nz,1,-1
           ic = (j-1)*nxz + (k-1)*nx + i
           IF (abs(fs_elev(ic)) > 0._kdp) THEN
              fs_elev1(ij) = fs_elev(ic)
           END IF
        END DO
     END DO
  END DO
  fs_elev = 0._kdp
  DO  i = 1,nx
     DO  j = 1,ny
        ij = (j-1)*nx + i
        IF(nz > 1) THEN
           DO k = ktop(ij),1,-1
              ic = (j-1)*nxz + (k-1)*nx + i        
              IF(p(ic) >= patm) THEN
                 k0 = k
                 k1 = k+1
                 IF(k0 < ktop(ij)) THEN
                    ! ... Bottom or intermediate plane
                    p0 = p(ic) - patm
                    p1 = p(ic+nx) - patm
                    z0 = zz(k0)
                    z1 = zz(k1)
                    fs_elev(ic) = (z1*p0-z0*p1)/(p0-p1)
                    EXIT
                 ELSEIF(k0 == ktop(ij)) THEN
                    ! ... Top plane
!                    p0 = p(ic-nx) - patm
                    p1 = p(ic) - patm
!                    z0 = zz(k0-1)
                    z1 = zz(k0)
                    fs_elev(ic) = p1/(denw(ic)*grav) + z1
                    EXIT
                 END IF
              END IF
           END DO
        ELSE
           ic = (j-1)*nxz + i        
           fs_elev(ic) = (p(ic)-patm)/(denw(ic)*grav) + zz(1)
        END IF
     END DO
  END DO
  ! ... Calculate the L2-norm of the change in water-table elevation
  DO  j = 1,ny
     DO  i = 1,nx
        ij = (j-1)*nx + i
        DO  k = nz,1,-1
           ic = (j-1)*nxz + (k-1)*nx + i
           IF (abs(fs_elev(ic)) > 0._kdp) THEN
              fs_elev2(ij) = fs_elev(ic)
           END IF
        END DO
     END DO
  END DO
  del_fs = fs_elev2 - fs_elev1
  dfs_l2norm = sqrt(dot_product(del_fs,del_fs))
  DEALLOCATE(fs_elev1, fs_elev2, del_fs)
END SUBROUTINE wt_elev
