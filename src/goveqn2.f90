SUBROUTINE goveqn2(gvm, gve, gveadv, gvecnd, nnodes, i, j, k)
  ! ... Purpose:  To sum the finite difference advective and conductive contributions to the
  ! ...      flow and heat transport governing equations for a given node
  ! ...      The returned values are split into advective and conductive parts for the
  ! ...      heat transport equation.
  USE machine_constants, ONLY: kdp
  USE bc, ONLY: ibc
  USE fdeq
  USE mesh
  USE parameters
  USE variables
  IMPLICIT NONE
  REAL(KIND=kdp), INTENT(OUT) :: gvm, gve, gveadv, gvecnd
  INTEGER, INTENT(OUT) :: nnodes
  INTEGER, INTENT(IN) :: i, j, k
  ! ...   GVM  - sum of mass governing eqn for node I,J,K
  ! ...   GVE  - sum of energy in governing eqn for node I,J,K
  ! ...   GVEADV - sum of advective energy governing eqn for node I,J,K
  ! ...   GVECND - sum of conductive energy governing eqn for node I,J,K
  ! ...   NNODES - number of active nodes adjacent to node IJK
  ! ...   I,J,K - indices of node
  !
  INTEGER :: ic, icxm, icxp, icym, icyp, iczm, iczp, ik, ikxm, ikxp,  &
       ikzm, ikzp, mxm, mxp, mym, myp, mzm, mzp
  REAL(KIND=kdp) :: ga, gc, ge, gm
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4 $//$Date: 2007/08/17 19:07:26 $'
  !     ------------------------------------------------------------------
  gvm = 0._kdp
  gve = 0._kdp
  gveadv = 0._kdp
  gvecnd = 0._kdp
  nnodes = 0
  ic = (k-1)*nx + i + (j-1)*nxz
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
        nnodes = nnodes + 1
        gvm = tx(mxp)*  &
             ((rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
             (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*(p(icxp)-p(ic))
        gveadv = tx(mxp)*  &
             ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
             (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))* &
             (p(icxp)-p(ic))
        gvecnd = txk(mxp)*(tc(icxp)-tc(ic))
        gve = gveadv + gvecnd
     END IF
  END IF
  ! ...           skip if this is the first node in the X-row
  ! ...           or is to the right of a non-active node
  IF (i > 1) THEN
     IF (ibc(ikxm,j) /= -1) THEN
        nnodes = nnodes + 1
        gm = tx(mxm)*((rs(ic)*usx(mxm)+rs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
             (rw(ic)*uwx(mxm)+rw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*(p(icxm)-p(ic))
        gvm = gvm + gm
        ga = tx(mxm)*  &
             ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
             (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*  &
             (p(icxm)-p(ic))
        gc = txk(mxm)*(tc(icxm)-tc(ic))
        ge = ga + gc
        gveadv = gveadv + ga
        gvecnd = gvecnd + gc
        gve = gve + ge
     END IF
  END IF
  ! ...           skip if this is the last node in the Y-row
  ! ...           or is to the left of a non-active node
  IF (j < ny) THEN
     IF (ibc(ik,j+1) /= -1) THEN
        nnodes = nnodes + 1
        gm = ty(myp)*((rs(icyp)*usy(myp)+rs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
             (rw(icyp)*uwy(myp)+rw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*(p(icyp)-p(ic))
        gvm = gvm + gm
        ga = ty(myp)*  &
             ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
             (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*  &
             (p(icyp)-p(ic))
        gc = tyk(myp)*(tc(icyp)-tc(ic))
        ge = ga + gc
        gveadv = gveadv + ga
        gvecnd = gvecnd + gc
        gve = gve + ge
     END IF
  END IF
  ! ...           skip if this is the first node in the Y-row
  ! ...           or is to the right of a non-active node
  IF (j > 1) THEN
     IF (ibc(ik,j-1) /= -1) THEN
        nnodes = nnodes + 1
        gm = ty(mym)*((rs(ic)*usy(mym)+rs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
             (rw(ic)*uwy(mym)+rw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*(p(icym)-p(ic))
        gvm = gvm + gm
        ga = ty(mym)*  &
             ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
             (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*  &
             (p(icym)-p(ic))
        gc = tyk(mym)*(tc(icym)-tc(ic))
        ge = ga + gc
        gveadv = gveadv + ga
        gvecnd = gvecnd + gc
        gve = gve + ge
     END IF
  END IF
  ! ...           skip if this is the top node in the Z-column
  ! ...           or is below a non-active node
  IF (k < nz) THEN
     IF (ibc(ikzp,j) /= -1) THEN
        nnodes = nnodes + 1
        gm = tz(mzp)*  &
             (((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
             (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
             (p(iczp)-p(ic))+  &
             ((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
             (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))
        gvm = gvm + gm
        ga = tz(mzp)*  &
             (((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
             (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
             (p(iczp)-p(ic))+  &
             ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
             (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))
        gc = tzk(mzp)*(tc(iczp)-tc(ic))
        ge = ga + gc
        gveadv = gveadv + ga
        gvecnd = gvecnd + gc
        gve = gve + ge
     END IF
  END IF
  ! ...           skip if this is the bottom node in the Z-column
  ! ...           or is above a non-active node
  IF (k > 1) THEN
     IF (ibc(ikzm,j) /= -1) THEN
        nnodes = nnodes + 1
        gm = tz(mzm)*  &
             (((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
             (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
             (p(iczm)-p(ic))-  &
             ((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
             (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))
        gvm = gvm + gm
        ga = tz(mzm)*  &
             (((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
             (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
             (p(iczm)-p(ic))-  &
             ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
             (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))
        gc = tzk(mzm)*(tc(iczm)-tc(ic))
        ge = ga + gc
        gveadv = gveadv + ga
        gvecnd = gvecnd + gc
        gve = gve + ge
     END IF
  END IF
END SUBROUTINE goveqn2
