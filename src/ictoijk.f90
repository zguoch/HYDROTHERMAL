SUBROUTINE ictoijk(ic,i,j,k,nx,nz)
  ! ... Returns the index (I,J,K) of the point with
  ! ...      vertical plane index IC.
  INTEGER, INTENT(IN) :: ic
  INTEGER, INTENT(OUT):: i, j, k
  INTEGER, INTENT(IN) :: nx, nz
  !
  INTEGER :: nxz
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5 $//$Date: 2006/03/07 23:45:45 $'
  !     ------------------------------------------------------------------
  !...
  nxz = nx*nz
  j = INT(REAL(ic-1)/REAL(nxz)) + 1
  k = INT(REAL((ic-1)-((j-1)*nxz))/REAL(nx)) + 1
  i = ic - ((k-1)*nx+(j-1)*nxz)
END SUBROUTINE ictoijk
