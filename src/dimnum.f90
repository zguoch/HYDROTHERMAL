SUBROUTINE dimnum
  ! ... Computes dimensionless numbers
  ! ...      Cell Thermal Peclet No.; Cell Nusselt No.
  ! ... Uses all available components of mass flux and temperature gradient
  USE machine_constants, ONLY: kdp
  USE bc, ONLY: ibc
  USE fdeq
  USE mesh
  USE parameters
  USE variables
  IMPLICIT NONE
  !
  INTEGER :: i, ic, icxm, icxp, icym, icyp, iczm, iczp, ik, ikxm, ikxp,  &
       ikzm, ikzp, j, k
  REAL(kind=kdp) :: alpham, dncpm, pex, pey, pez,  &
       ptpx, ptpxm, ptpxp, ptpy, ptpym, ptpyp, ptpz, ptpzm, ptpzp, wt,  &
       qeadv, qecond
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.2 $//$Date: 2002/10/28 16:26:25 $'
  !     ------------------------------------------------------------------
  ! ... Cell Nusselt Number
  nu = 0._kdp
  DO ic=1,nxyz
     qeadv = 0._kdp
     qecond = 0._kdp
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     ! ... Test for nodes and neighbors in active region
     !              X direction
     IF (i > 1 .AND. i < nx) THEN
        icxp = ic + 1
        icxm = ic - 1
        ikxp = ik + 1
        ikxm = ik - 1
        IF (ibc(ik,j) /= -1 .AND. ibc(ikxp,j) /= -1 .AND. ibc(ikxm,j) /= -1) THEN
           ptpxp = (tc(icxp) - tc(ic))/(xx(i+1) - xx(i))
           ptpxm = (tc(ic) - tc(icxm))/(xx(i) - xx(i-1))
           wt = (xx(i) - xx(i-1))/(xx(i+1) - xx(i-1))
           ptpx = wt*ptpxp + (1._kdp-wt)*ptpxm
        END IF
     END IF
     !              Y direction
     IF (j > 1 .AND. j < ny) THEN
        icyp = ic + nxz
        icym = ic - nxz
        IF (ibc(ik,j) /= -1 .AND. ibc(ik,j+1) /= -1 .AND. ibc(ik,j-1) /= -1) THEN
           ptpyp = (tc(icyp) - tc(ic))/(yy(j+1) - yy(j))
           ptpym = (tc(ic) - tc(icym))/(yy(j) - yy(j-1))
           wt = (yy(j) - yy(j-1))/(yy(j+1) - yy(j-1))
           ptpy = wt*ptpyp + (1._kdp-wt)*ptpym
        END IF
     END IF
     !              Z direction
     IF (k > 1 .AND. k < nz) THEN
        iczp = ic + nx
        iczm = ic - nx
        ikzp = ik + nx
        ikzm = ik - nx
        IF (ibc(ik,j) /= -1 .AND. ibc(ikzp,j) /= -1 .AND. ibc(ikzm,j) /= -1) THEN
           ptpzp = (tc(iczp) - tc(ic))/(zz(k+1) - zz(k))
           ptpzm = (tc(ic) - tc(iczm))/(zz(k) - zz(k-1))
           wt = (zz(k) - zz(k-1))/(zz(k+1) - zz(k-1))
           ptpz = wt*ptpzp + (1._kdp-wt)*ptpzm
        END IF
     END IF
     qecond = xkc(ic)*sqrt(ptpx*ptpx+ptpy*ptpy+ptpz*ptpz)   ! ... Only isotropic conduction
     IF(qecond <= 0._kdp) CYCLE
     qeadv = h(ic)*sqrt(xwmflx(ic)*xwmflx(ic)+ywmflx(ic)*ywmflx(ic)+zwmflx(ic)*zwmflx(ic)+  &
          xsmflx(ic)*xsmflx(ic)+ysmflx(ic)*ysmflx(ic)+zsmflx(ic)*zsmflx(ic))
     nu(ic) = (qeadv+qecond)/qecond
  END DO
  ! ... Cell Peclet Number
  ! ... Only for single phase water cells at present 
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     ! ... skip 2 phase cells and air-water cells and steam cells
     IF(ibc(ik,j) /= -1 .AND. ind(ic) /=2 .AND. ind(ic) /= 3 .AND. ind(ic) /= 5) THEN
        ! ... Cell Peclet number as rms of directional Peclet numbers
        dncpm = phi(ic)*denw(ic)/dth(ic) + (1._kdp-phi(ic))*df(ic)*phfwt(ic)
        alpham = xkc(ic)/dncpm     ! ... Only isotropic medium conductivity available
        pex = (dx(i)*(xwmflx(ic)/denw(ic)))/alpham
        pey = (dy(j)*(ywmflx(ic)/denw(ic)))/alpham
        pez = (dz(k)*(zwmflx(ic)/denw(ic)))/alpham
        pe(ic) = sqrt(pex*pex+pey*pey+pez*pez)
     END IF
  END DO
END SUBROUTINE dimnum
