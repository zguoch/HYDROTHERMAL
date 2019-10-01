SUBROUTINE goveqn(gveqnm, gveqne, i, j, k)
  ! ... Purpose:  To sum the finite difference advective and conductive contributions to the
  ! ...        flow and heat transport governing equations for a given node
  USE machine_constants, ONLY: kdp
  USE bc, ONLY: ibc
  USE fdeq
  USE mesh
  USE parameters
  USE variables
  IMPLICIT NONE
  REAL(KIND=kdp), INTENT(OUT) :: gveqnm, gveqne
  INTEGER, INTENT(IN) :: i, j, k
  ! ...   GVEQNM - sum of mass governing eqn for node I,J,K
  ! ...   GVEQNE - sum of energy governing eqn for node I,J,K
  ! ...   I,J,K - indices of node
  !
  INTEGER :: ic, icxm, icxp, icym, icyp, iczm, iczp, ik, ikxm, ikxp,  &
       ikzm, ikzp, mxm, mxp, mym, myp, mzm
  INTEGER :: mzp
  REAL(KIND=kdp) :: gove, govm
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4 $//$Date: 2007/08/24 21:40:28 $'
  !     ------------------------------------------------------------------
  !...
  gveqnm = 0._kdp
  gveqne = 0._kdp
  govm = 0._kdp
  gove = 0._kdp
  ic = (k-1)*nx + i + (j-1)*nx*nz
  ik = (k-1)*nx + i
  ! ...          X direction
  icxp = ic + 1
  icxm = ic - 1
  ikxp = ik + 1
  ikxm = ik - 1
  ! ...          Z direction
  iczp = ic + nx
  iczm = ic - nx
  ikzp = ik + nx
  ikzm = ik - nx
  ! ...          Y direction
  icyp = ic + nx*nz
  icym = ic - nx*nz
  ! ...       intercell regions X direction
  mxm = (k-1)*nxx + i + (j-1)*nxx*nz
  mxp = mxm + 1
  ! ...       intercell regions Z direction
  mzm = (k-1)*nx + i + (j-1)*nx*nzz
  mzp = mzm + nx
  ! ...       intercell regions Y direction
  mym = ic
  myp = mym + nx*nz
  ! ... Positive values for gveqn_ mean inflow to the cell
  ! ...           skip if this is the last node in the X-row
  ! ...           or is to the left of a non-active node
  IF (i < nx) THEN
     IF (ibc(ikxp,j) /= -1) THEN
        gveqnm = gveqnm + tx(mxp)*  &
             ((rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
             (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*(p(icxp)-p(ic))
        gveqne = gveqne + tx(mxp)*  &
             ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
             (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*  &
             (p(icxp)-p(ic)) + txk(mxp)*(tc(icxp)-tc(ic))
     END IF
  END IF
  ! ...           skip if this is the first node in the X-row
  ! ...           or is to the right of a non-active node
  IF (i > 1) THEN
     IF (ibc(ikxm,j) /= -1) THEN
        govm = tx(mxm)*  &
             ((rs(ic)*usx(mxm)+rs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
             (rw(ic)*uwx(mxm)+rw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*(p(icxm)-p(ic))
        gveqnm = gveqnm + govm
        gove = tx(mxm)*  &
             ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
             (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*  &
             (p(icxm)-p(ic)) + txk(mxm)*(tc(icxm)-tc(ic))
        gveqne = gveqne + gove
     END IF
  END IF
  ! ...           skip if this is the last node in the Y-row
  ! ...           or is to the left of a non-active node
  IF (j < ny) THEN
     IF (ibc(ik,j+1) /= -1) THEN
        govm = ty(myp)*  &
             ((rs(icyp)*usy(myp)+rs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
             (rw(icyp)*uwy(myp)+rw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*(p(icyp)-p(ic))
        gveqnm = gveqnm + govm
        gove = ty(myp)*  &
             ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
             (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*  &
             (p(icyp)-p(ic)) + tyk(myp)*(tc(icyp)-tc(ic))
        gveqne = gveqne + gove
     END IF
  END IF
  ! ...           skip if this is the first node in the Y-row
  ! ...           or is to the right of a non-active node
  IF (j > 1) THEN
     IF (ibc(ik,j-1) /= -1) THEN
        govm = ty(mym)*  &
             ((rs(ic)*usy(mym)+rs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
             (rw(ic)*uwy(mym)+rw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*(p(icym)-p(ic))
        gveqnm = gveqnm + govm
        gove = ty(mym)*  &
             ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
             (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*  &
             (p(icym)-p(ic)) + tyk(mym)*(tc(icym)-tc(ic))
        gveqne = gveqne + gove
     END IF
  END IF
  ! ...           skip if this is the top node in the Z-column
  ! ...           or is below a non-active node
  IF (k < nz) THEN
     IF (ibc(ikzp,j) /= -1) THEN
        govm = tz(mzp)*  &
             (((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
             (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
             (p(iczp)-p(ic))+  &
             ((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
             (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))
        gveqnm = gveqnm + govm
        gove = tz(mzp)*  &
             (((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
             (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
             (p(iczp)-p(ic))+  &
             ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
             (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))  &
             + tzk(mzp)*(tc(iczp)-tc(ic))
        gveqne = gveqne + gove
     END IF
  END IF
  ! ...           skip if this is the bottom node in the Z-column
  ! ...           or is above a non-active node
  IF (k > 1) THEN
     IF (ibc(ikzm,j) /= -1) THEN
        govm = tz(mzm)*  &
             (((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
             (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
             (p(iczm)-p(ic))-  &
             ((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
             (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))
        gveqnm = gveqnm + govm
        gove = tz(mzm)*  &
             (((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
             (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
             (p(iczm)-p(ic))-  &
             ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
             (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))  &
             + tzk(mzm)*(tc(iczm)-tc(ic))
        gveqne = gveqne + gove
     END IF
  END IF
END SUBROUTINE goveqn
