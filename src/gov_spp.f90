SUBROUTINE gov_spp(gveqnm, i, j, k)
  ! ... Purpose:  To sum the finite difference conductive flow contributions to the
  ! ...        governing equation of flow for a given node.
  ! ... For specified pressure nodes to compute inflow or outflow rate
  USE machine_constants, ONLY: kdp
!!$  USE bc, ONLY: ibc
  USE fdeq
  USE mesh
  USE parameters
  USE variables
  IMPLICIT NONE
  REAL(KIND=kdp), INTENT(OUT) :: gveqnm
  INTEGER, INTENT(IN) :: i, j, k
  ! ...   GVEQNM - sum of mass governing eqn for node I,J,K
  ! ...   I,J,K - indices of node
  !
  INTEGER :: ic, icxm, icxp, icym, icyp, iczm, iczp, ik, ikxm, ikxp,  &
       ikzm, ikzp, mxm, mxp, mym, myp, mzm, mzp
  REAL(KIND=kdp) :: govm
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3 $//$Date: 2007/08/17 19:07:26 $'
  !     ------------------------------------------------------------------
  !...
  ic = (k-1)*nx + i + (j-1)*nxz
  ik = (k-1)*nx + i
  icxp = ic + 1
  icxm = ic - 1
  ikxp = ik + 1
  ikxm = ik - 1
  iczp = ic + nx
  iczm = ic - nx
  ikzp = ik + nx
  ikzm = ik - nx
  icyp = ic + nxz
  icym = ic - nxz
  ! ... cell face indices
  mxm = (k-1)*nxx + i + (j-1)*nxx*nz
  mxp = mxm + 1
  mzm = (k-1)*nx + i + (j-1)*nx*nzz
  mzp = mzm + nx
  mym = ic
  myp = mym + nxz
  gveqnm = 0._kdp
  govm = 0._kdp
  ! ... Positive value for gveqnm means inflow to the cell
  ! ...           skip if this is the last node in the X-row
  ! ...      depends on t_ being zero for boundary cell faces
  IF (i < nx) THEN
     gveqnm = gveqnm + tx(mxp)*  &
          ((rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
          (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*(p(icxp)-p(ic))
  END IF
  ! ...           skip if this is the first node in the X-row
  IF (i > 1) THEN
     govm = tx(mxm)*  &
          ((rs(ic)*usx(mxm)+rs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
          (rw(ic)*uwx(mxm)+rw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*(p(icxm)-p(ic))
     gveqnm = gveqnm + govm
  END IF
  ! ...           skip if this is the last node in the Y-row
  IF (j < ny) THEN
     govm = ty(myp)*  &
          ((rs(icyp)*usy(myp)+rs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
          (rw(icyp)*uwy(myp)+rw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*(p(icyp)-p(ic))
     gveqnm = gveqnm + govm
  END IF
  ! ...           skip if this is the first node in the Y-row
  IF (j > 1) THEN
     govm = ty(mym)*  &
          ((rs(ic)*usy(mym)+rs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
          (rw(ic)*uwy(mym)+rw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*(p(icym)-p(ic))
     gveqnm = gveqnm + govm
  END IF
  ! ...           skip if this is the top node in the Z-column
  IF (k < nz) THEN
     govm = tz(mzp)*  &
          (((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
          (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
          (p(iczp)-p(ic))+  &
          ((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
          (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))
     gveqnm = gveqnm + govm
  END IF
  ! ...           skip if this is the bottom node in the Z-column
  IF (k > 1) THEN
     govm = tz(mzm)*  &
          (((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
          (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
          (p(iczm)-p(ic))-  &
          ((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
          (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))
     gveqnm = gveqnm + govm
  END IF
END SUBROUTINE gov_spp
