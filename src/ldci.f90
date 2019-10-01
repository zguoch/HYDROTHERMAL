SUBROUTINE ldci
  ! ... Loads the CI array for connected nodes based on the renumbered
  ! ...      mesh of active nodes
  ! ... The CI array contains the node numbers of the active nodes connected
  ! ...      to a given active node
  USE mesh
  INTEGER :: ic, i, j, k, jc, m, ma
  INTEGER, DIMENSION(6) :: icc
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5 $//$Date: 2007/08/08 17:06:19 $'
  !     ------------------------------------------------------------------
  DO  m=1,nxyz
     CALL mtoijk(m,i,j,k,nx,ny)
     ic = (j-1)*nx*nz + (k-1)*nx + i
     ma = mrno(ic)
     ! ... GMRES solver: MA is non zero for all active cells in the rectangular mesh
     ! ... SSOR  solver: MA is non zero for all active cells not specified value b.c.
     ! ...      in the rectangular mesh
     ! ...      Ordered by x-z planes
     IF(ma > 0) THEN
        icc(1) = ic-nx
        icc(2) = ic-nx*nz
        icc(3) = ic-1
        icc(4) = ic+1
        icc(5) = ic+nx*nz
        icc(6) = ic+nx
        ! ... ICC is a full global mesh node number in x,z,y order
        ! ...     Some values will be out of range for nodes on boundaries
        DO  jc=1,6
           IF ((icc(jc) >= 1) .AND. (icc(jc) <= nxyz)) ci(jc,ma) = mrno(icc(jc))
        END DO
        ! ... Modify indices for ends of node rows; fix wrap around effect
        IF (k == 1) ci(1,ma) = 0
        IF (j == 1) ci(2,ma) = 0
        IF (i == 1) ci(3,ma) = 0
        IF (i == nx) ci(4,ma) = 0
        IF (j == ny) ci(5,ma) = 0
        IF (k == nz) ci(6,ma) = 0
     END IF
  END DO
END SUBROUTINE ldci
