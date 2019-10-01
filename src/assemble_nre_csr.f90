SUBROUTINE assemble_nre_csr
  !     Purpose: To assemble the set of Newton-Raphson matrix equations.
  !              Based on the fully implicit finite difference
  !              equations for flow and energy transport and the four
  !              partial derivative equations (partials of the residual flow 
  !              and transport equations with respect to pressure and enthalpy).
  USE machine_constants, ONLY: kdp
  USE math_constants
  USE bc
  USE fdeq
  USE mesh
  USE parameters
  USE solver_gmres
  USE source
  USE variables
  IMPLICIT NONE
  ! 
  INTEGER :: i, ic, icxm, icxp, icym, icyp, iczm, iczp, ik, j, k, ma
  INTEGER :: mxm, mxp, mym, myp, mzm, mzp
  INTEGER :: nnn
  REAL(KIND=kdp) :: vol
  REAL(KIND=kdp) :: gveqnm, gveqne, qsppbc, qhspp, dqhhspp, dqhpspp,  &
       qseepbc, qhseep, dqhhseep, dqhpseep 
  REAL(KIND=kdp) :: ablkpp, ablkph, ablkhp, ablkhh,  &
       ablkpp1, ablkph1, ablkhp1, ablkhh1,  &
       ablkpp2, ablkph2, ablkhp2, ablkhh2,  &
       ablkpp3, ablkph3, ablkhp3, ablkhh3,  &
       ablkpp4, ablkph4, ablkhp4, ablkhh4,  &
       ablkpp5, ablkph5, ablkhp5, ablkhh5,  &
       ablkpp6, ablkph6, ablkhp6, ablkhh6,  &
       rhsp, rhsh
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.5 $//$Date: 2007/08/24 00:01:59 $'
  !     ------------------------------------------------------------------
  !...
  ! ... Uses compressed row storage of the sparse matrix with off-diagonal elements
  ! ...      in any order. Here the order is by local face number, [1-6].
  ! ... Main diagonal stored first in each row (not required by this gmres)
  ! ... Initialize arrays RHS and AV
  av = 0._kdp
  rhs = 0._kdp
  nnn = 0          ! ... initialize counter for array element
  ! ... Loop over all mesh nodes; Assemble for each active node into av and rhs arrays
  DO ic=1,nxyz
     ma = mrno(ic)
     IF(ma == 0) CYCLE
     ! ... Initialize for current block 2x2 sub-array
     ablkpp = 0._kdp
     ablkph = 0._kdp
     ablkhp = 0._kdp
     ablkhh = 0._kdp
     rhsp = 0._kdp
     rhsh = 0._kdp
     ! ... Adjacent node indices
     icxp = ic + 1
     icxm = ic - 1
     icyp = ic + nxz
     icym = ic - nxz
     iczp = ic + nx
     iczm = ic - nx
     CALL ictoijk(ic,i,j,k,nx,nz)
     ! ... Node number in current slice, j
     ik = (k-1)*nx + i
     ! ... Cell face indices for conductance factors
     mxm = (k-1)*nxx + i + (j-1)*nxx*nz
     mxp = mxm + 1
     mym = ic
     myp = mym + nxz
     mzm = (k-1)*nx + i + (j-1)*nx*nzz
     mzp = mzm + nx
     IF(ibc(ik,j) == 11) THEN
        ! ... Assembly of specified pressure cell with specified enthalpy
        !    ******* Main Diagonal Block Terms and Right Hand Side *********
        ! ... Flow N-R equation
        ! ...      Assemble a trivial equation for pressure 
        ! ... Heat N-R equation
        ! ...      Assemble a trivial equation for enthalpy
        ! ... main diagonal stored first in each row
        ablkpp = -1._kdp
        ablkph = 0._kdp
        ablkhp = 0._kdp
        ablkhh = -1._kdp
        rhsp = 0._kdp
        rhsh = 0._kdp
        !    ******* Off Diagonal Block Terms  ***************
        ! ... derivatives of residual flow and transport equations w.r.t. P & H at l+1
        IF(ci(4,ma) > 0) THEN
           ablkpp4 = 0._kdp
           ablkph4 = 0._kdp
           ablkhp4 = 0._kdp
           ablkhh4 = 0._kdp
        END IF
        IF(ci(3,ma) > 0) THEN
           ablkpp3 = 0._kdp
           ablkph3 = 0._kdp
           ablkhp3 = 0._kdp
           ablkhh3 = 0._kdp
        END IF
        IF(ci(5,ma) > 0) THEN
           ablkpp5 = 0._kdp
           ablkph5 = 0._kdp
           ablkhp5 = 0._kdp
           ablkhh5 = 0._kdp
        END IF
        IF(ci(2,ma) > 0) THEN
           ablkpp2 = 0._kdp
           ablkph2 = 0._kdp
           ablkhp2 = 0._kdp
           ablkhh2 = 0._kdp
        END IF
        IF(ci(6,ma) > 0) THEN
           ablkpp6 = 0._kdp
           ablkph6 = 0._kdp
           ablkhp6 = 0._kdp
           ablkhh6 = 0._kdp
        END IF
        IF(ci(1,ma) > 0) THEN
           ablkpp1 = 0._kdp
           ablkph1 = 0._kdp
           ablkhp1 = 0._kdp
           ablkhh1 = 0._kdp
        END IF
     ELSEIF(ibc(ik,j) == 10) THEN
        ! ... Assembly of specified pressure cell with associated enthalpy
        !    ******* Main Diagonal Block Terms and Right Hand Side *********
        ! ... Flow N-R equation
        ! ...      Assemble a trivial equation for pressure 
        ablkpp = -1._kdp
        ablkph = 0._kdp
        rhsp = 0._kdp
        ! ... Heat
        IF(ci(4,ma) > 0) THEN
           ablkhp = ablkhp - tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                (dxvdsp(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsp(mxp)*hrs(ic)+dhrsp(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwp(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwp(mxp)*hrw(ic)+dhrwp(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp))) - txk(mxp)*dtp(ic)
           ablkhh = ablkhh - tx(mxp)*(dxvdsh(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsh(mxp)*hrs(ic)+dhrsh(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwh(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwh(mxp)*hrw(ic)+dhrwh(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp)) - txk(mxp)*dth(ic)
           rhsh = rhsh - tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*  &
                (p(icxp)-p(ic)) - txk(mxp)*(tc(icxp)-tc(ic))
        END IF
        IF(ci(3,ma) > 0) THEN
           ablkhp = ablkhp - tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                ((dxvdsp(mxm)*hrs(ic)+dhrsp(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsp(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwp(mxm)*hrw(ic)+dhrwp(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwp(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm)))  &
                - txk(mxm)*dtp(ic)
           ablkhh = ablkhh - tx(mxm)*((dxvdsh(mxm)*hrs(ic)+dhrsh(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsh(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwh(mxm)*hrw(ic)+dhrwh(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwh(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm))  &
                - txk(mxm)*dth(ic)
           rhsh = rhsh - tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*  &
                (p(icxm)-p(ic)) - txk(mxm)*(tc(icxm)-tc(ic))
        END IF
        IF(ci(5,ma) > 0) THEN
           ablkhp = ablkhp - ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                (dyvdsp(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsp(myp)*hrs(ic)+dhrsp(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwp(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwp(myp)*hrw(ic)+dhrwp(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp))) - tyk(myp)*dtp(ic)
           ablkhh = ablkhh - ty(myp)*(dyvdsh(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsh(myp)*hrs(ic)+dhrsh(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwh(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwh(myp)*hrw(ic)+dhrwh(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp)) - tyk(myp)*dth(ic)
           rhsh = rhsh - ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*  &
                (p(icyp)-p(ic)) - tyk(myp)*(tc(icyp)-tc(ic))
        END IF
        IF(ci(2,ma) > 0) THEN
           ablkhp = ablkhp - ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                ((dyvdsp(mym)*hrs(ic)+dhrsp(ic)*yvds(mym))*usy(mym)+  &
                dyvdsp(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwp(mym)*hrw(ic)+dhrwp(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwp(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym)))  &
                - tyk(mym)*dtp(ic)
           ablkhh = ablkhh - ty(mym)*((dyvdsh(mym)*hrs(ic)+dhrsh(ic)*yvds(mym))*usy(mym)+  &
                dyvdsh(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwh(mym)*hrw(ic)+dhrwh(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwh(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym))  &
                - tyk(mym)*dth(ic)
           rhsh = rhsh -ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*  &
                (p(icym)-p(ic)) - tyk(mym)*(tc(icym)-tc(ic))
        END IF
        IF(ci(6,ma) > 0) THEN
           ablkhp = ablkhp - tz(mzp)*  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp)+  &
                (dzvdsp(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdsp(mzp)*hrs(ic)+dhrsp(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwp(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdwp(mzp)*hrw(ic)+dhrwp(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsp(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdgsp(mzp)*hrs(ic)+dhrsp(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwp(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdgwp(mzp)*hrw(ic)+dhrwp(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k)) - tzk(mzp)*dtp(ic)
           ablkhh = ablkhh - tz(mzp)*((dzvdsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdsh(mzp)*hrs(ic)+dhrsh(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdwh(mzp)*hrw(ic)+dhrwh(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdgsh(mzp)*hrs(ic)+dhrsh(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdgwh(mzp)*hrw(ic)+dhrwh(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k)) - tzk(mzp)*dth(ic)
           rhsh = rhsh - tz(mzp)*  &
                (((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
                (p(iczp)-p(ic))+  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))  &
                - tzk(mzp)*(tc(iczp)-tc(ic))
        END IF
        IF(ci(1,ma) > 0) THEN
           ablkhp = ablkhp - tz(mzm)*  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm)+  &
                ((dzvdsp(mzm)*hrs(ic)+dhrsp(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsp(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwp(mzm)*hrw(ic)+dhrwp(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwp(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsp(mzm)*hrs(ic)+dhrsp(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsp(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwp(mzm)*hrw(ic)+dhrwp(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwp(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))  &
                - tzk(mzm)*dtp(ic)
           ablkhh = ablkhh - tz(mzm)*(((dzvdsh(mzm)*hrs(ic)+dhrsh(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwh(mzm)*hrw(ic)+dhrwh(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsh(mzm)*hrs(ic)+dhrsh(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwh(mzm)*hrw(ic)+dhrwh(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))  &
                - tzk(mzm)*dth(ic)
           rhsh = rhsh - tz(mzm)*  &
                (((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
                (p(iczm)-p(ic))-  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))  &
                - tzk(mzm)*(tc(iczm)-tc(ic))
        END IF
        ! ... Get specified pressure cell boundary face flow rate from residual and calculate 
        ! ...     advective heat gain or loss term for heat equation. 
        ! ...     It is lagged one NR iteration.
        CALL gov_spp(gveqnm, i, j, k)
        ! ... Add in the capacitance terms
        ! ...     Fluid source in specified pressure cell is not realistic
        gveqnm = gveqnm - (xm(ic)-xmoldt(ic))/delt
        qsppbc = -gveqnm     ! ... Residual net inflow is specified p cell outflow
        ! ... Calculate the heat b.c. source term for a specified p node
        ! ... Only the case of single-phase seepage is allowed.
        ! ... No other heat sources are allowed in the specified pressure cell
        IF(qsppbc > 0._kdp) THEN     ! ... Inflow
           qhspp = qsppbc*ehassoc(ic)
           dqhhspp = 0._kdp         ! ... Jacobian terms for heat eq.
        ELSE                         ! ... Outflow
           qhspp = qsppbc*h(ic)
           dqhhspp = qsppbc         ! ... Jacobian terms for heat eq.
        END IF
        dqhpspp = 0._kdp     ! ... Constant pressure at cell so dqhp*delp=0
        vol = dx(i)*dy(j)*dz(k)
        IF(irad) vol = (rsq(i+1)-rsq(i))*pi*dz(k)
        ablkhp = ablkhp - c(ic)*vol/delt + dqhpspp
        ablkhh = ablkhh - d(ic)*vol/delt + dqhhspp
        ! ... Add in the capacitance and b.c.source terms
        rhsh = rhsh + (en(ic)-enoldt(ic))/delt - qhspp
        !    ******* Off Diagonal Block Terms  ***************
        ! ... derivatives of residual flow and transport equations w.r.t. P & H at l+1
        IF(ci(4,ma) > 0) THEN
           ablkpp4 = 0._kdp
           ablkph4 = 0._kdp
           ablkhp4 = tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                ((dxvdsp(mxp)*hrs(icxp)+dhrsp(icxp)*xvds(mxp))*usx(mxp)+  &
                dxvdsp(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                (dxvdwp(mxp)*hrw(icxp)+dhrwp(icxp)*xvdw(mxp))*uwx(mxp)+  &
                dxvdwp(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic)))  &
                + txk(mxp)*dtp(icxp)
           ablkhh4 = tx(mxp)*  &
                ((dxvdsh(mxp)*hrs(icxp)+dhrsh(icxp)*xvds(mxp))*usx(mxp)+  &
                dxvdsh(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                (dxvdwh(mxp)*hrw(icxp)+dhrwh(icxp)*xvdw(mxp))*uwx(mxp)+  &
                dxvdwh(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic))  &
                + txk(mxp)*dth(icxp)
        END IF
        IF(ci(3,ma) > 0) THEN
           ablkpp3 = 0._kdp
           ablkph3 = 0._kdp
           ablkhp3 = tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                (dxvdsp(mxm)*hrs(ic)*usx(mxm)+  &
                (dxvdsp(mxm)*hrs(icxm)+dhrsp(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                dxvdwp(mxm)*hrw(ic)*uwx(mxm)+  &
                (dxvdwp(mxm)*hrw(icxm)+dhrwp(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                (p(icxm)-p(ic))) + txk(mxm)*dtp(icxm)
           ablkhh3 = tx(mxm)*(dxvdsh(mxm)*hrs(ic)*usx(mxm)+  &
                (dxvdsh(mxm)*hrs(icxm)+dhrsh(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                dxvdwh(mxm)*hrw(ic)*uwx(mxm)+  &
                (dxvdwh(mxm)*hrw(icxm)+dhrwh(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                (p(icxm)-p(ic)) + txk(mxm)*dth(icxm)
        END IF
        IF(ci(5,ma) > 0) THEN
           ablkpp5 = 0._kdp
           ablkph5 = 0._kdp
           ablkhp5 = ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                ((dyvdsp(myp)*hrs(icyp)+dhrsp(icyp)*yvds(myp))*usy(myp)+  &
                dyvdsp(myp)*hrs(ic)*(1._kdp-usy(myp))+  &
                (dyvdwp(myp)*hrw(icyp)+dhrwp(icyp)*yvdw(myp))*uwy(myp)+  &
                dyvdwp(myp)*hrw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic)))  &
                + tyk(myp)*dtp(icyp)
           ablkhh5 = ty(myp)*  &
                ((dyvdsh(myp)*hrs(icyp)+dhrsh(icyp)*yvds(myp))*usy(myp)+  &
                dyvdsh(myp)*hrs(ic)*(1._kdp-usy(myp))+  &
                (dyvdwh(myp)*hrw(icyp)+dhrwh(icyp)*yvdw(myp))*uwy(myp)+  &
                dyvdwh(myp)*hrw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic))  &
                + tyk(myp)*dth(icyp)
        END IF
        IF(ci(2,ma) > 0) THEN
           ablkpp2 = 0._kdp
           ablkph2 = 0._kdp
           ablkhp2 = ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                (dyvdsp(mym)*hrs(ic)*usy(mym)+  &
                (dyvdsp(mym)*hrs(icym)+dhrsp(icym)*yvds(mym))*(1._kdp-usy(mym))+  &
                dyvdwp(mym)*hrw(ic)*uwy(mym)+  &
                (dyvdwp(mym)*hrw(icym)+dhrwp(icym)*yvdw(mym))*(1._kdp-uwy(mym)))*  &
                (p(icym)-p(ic))) + tyk(mym)*dtp(icym)
           ablkhh2 = ty(mym)*(dyvdsh(mym)*hrs(ic)*usy(mym)+  &
                (dyvdsh(mym)*hrs(icym)+dhrsh(icym)*yvds(mym))*(1._kdp-usy(mym))+  &
                dyvdwh(mym)*hrw(ic)*uwy(mym)+  &
                (dyvdwh(mym)*hrw(icym)+dhrwh(icym)*yvdw(mym))*(1._kdp-uwy(mym)))*  &
                (p(icym)-p(ic)) + tyk(mym)*dth(icym)
        END IF
        IF(ci(6,ma) > 0) THEN
           ablkpp6 = 0._kdp
           ablkph6 = 0._kdp
           ablkhp6 = tz(mzp)*  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp)+  &
                ((dzvdsp(mzp)*hrs(iczp)+dhrsp(iczp)*zvds(mzp))*usz(mzp)+  &
                dzvdsp(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdwp(mzp)*hrw(iczp)+dhrwp(iczp)*zvdw(mzp))*uwz(mzp)+  &
                dzvdwp(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                ((dzvdgsp(mzp)*hrs(iczp)+dhrsp(iczp)*zvdgs(mzp))*usz(mzp)+  &
                dzvdgsp(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdgwp(mzp)*hrw(iczp)+dhrwp(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                dzvdgwp(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*zzz(k))  &
                + tzk(mzp)*dtp(iczp)
           ablkhh6 = tz(mzp)*  &
                (((dzvdsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvds(mzp))*usz(mzp)+  &
                dzvdsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdw(mzp))*uwz(mzp)+  &
                dzvdwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                ((dzvdgsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvdgs(mzp))*usz(mzp)+  &
                dzvdgsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdgwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                dzvdgwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*zzz(k))  &
                + tzk(mzp)*dth(iczp)
        END IF
        IF(ci(1,ma) > 0) THEN
           ablkpp1 = 0._kdp
           ablkph1 = 0._kdp
           ablkhp1 = tz(mzm)*  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm)+  &
                (dzvdsp(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdsp(mzm)*hrs(iczm)+dhrsp(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                dzvdwp(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdwp(mzm)*hrw(iczm)+dhrwp(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                (p(iczm)-p(ic))- (dzvdgsp(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdgsp(mzm)*hrs(iczm)+dhrsp(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                dzvdgwp(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdgwp(mzm)*hrw(iczm)+dhrwp(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                zzz(k-1)) + tzk(mzm)*dtp(iczm)
           ablkhh1 = tz(mzm)*((dzvdsh(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                dzvdwh(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                (p(iczm)-p(ic))- (dzvdgsh(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdgsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                dzvdgwh(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdgwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                zzz(k-1)) + tzk(mzm)*dth(iczm)
        END IF
     ELSEIF((ibc(ik,j)/10 == 3 .OR. ibc(ik,j)/10 == 5) .AND. seeping(ik,j)) THEN
        ! ... Assembly of seeping cell
        !    ******* Main Diagonal Block Terms and Right Hand Side *********
        ! ... Flow N-R equation
        ! ...      Assemble a trivial equation for pressure 
        ! ... main diagonal stored first in each row
        ablkpp = -1._kdp
        ablkph = 0._kdp
        rhsp = 0._kdp
        ! ... Heat
        IF(ci(4,ma) > 0) THEN
           ablkhp = ablkhp - tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                (dxvdsp(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsp(mxp)*hrs(ic)+dhrsp(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwp(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwp(mxp)*hrw(ic)+dhrwp(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp))) - txk(mxp)*dtp(ic)
           ablkhh = ablkhh - tx(mxp)*(dxvdsh(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsh(mxp)*hrs(ic)+dhrsh(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwh(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwh(mxp)*hrw(ic)+dhrwh(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp)) - txk(mxp)*dth(ic)
           rhsh = rhsh - tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*  &
                (p(icxp)-p(ic)) - txk(mxp)*(tc(icxp)-tc(ic))
        END IF
        IF(ci(3,ma) > 0) THEN
           ablkhp = ablkhp - tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                ((dxvdsp(mxm)*hrs(ic)+dhrsp(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsp(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwp(mxm)*hrw(ic)+dhrwp(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwp(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm)))  &
                - txk(mxm)*dtp(ic)
           ablkhh = ablkhh - tx(mxm)*((dxvdsh(mxm)*hrs(ic)+dhrsh(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsh(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwh(mxm)*hrw(ic)+dhrwh(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwh(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm))  &
                - txk(mxm)*dth(ic)
           rhsh = rhsh - tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*  &
                (p(icxm)-p(ic)) - txk(mxm)*(tc(icxm)-tc(ic))
        END IF
        IF(ci(5,ma) > 0) THEN
           ablkhp = ablkhp - ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                (dyvdsp(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsp(myp)*hrs(ic)+dhrsp(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwp(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwp(myp)*hrw(ic)+dhrwp(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp))) - tyk(myp)*dtp(ic)
           ablkhh = ablkhh - ty(myp)*(dyvdsh(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsh(myp)*hrs(ic)+dhrsh(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwh(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwh(myp)*hrw(ic)+dhrwh(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp)) - tyk(myp)*dth(ic)
           rhsh = rhsh - ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*  &
                (p(icyp)-p(ic)) - tyk(myp)*(tc(icyp)-tc(ic))
        END IF
        IF(ci(2,ma) > 0) THEN
           ablkhp = ablkhp - ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                ((dyvdsp(mym)*hrs(ic)+dhrsp(ic)*yvds(mym))*usy(mym)+  &
                dyvdsp(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwp(mym)*hrw(ic)+dhrwp(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwp(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym)))  &
                - tyk(mym)*dtp(ic)
           ablkhh = ablkhh - ty(mym)*((dyvdsh(mym)*hrs(ic)+dhrsh(ic)*yvds(mym))*usy(mym)+  &
                dyvdsh(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwh(mym)*hrw(ic)+dhrwh(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwh(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym))  &
                - tyk(mym)*dth(ic)
           rhsh = rhsh -ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*  &
                (p(icym)-p(ic)) - tyk(mym)*(tc(icym)-tc(ic))
        END IF
        IF(ci(6,ma) > 0) THEN
           ablkhp = ablkhp - tz(mzp)*  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp)+  &
                (dzvdsp(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdsp(mzp)*hrs(ic)+dhrsp(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwp(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdwp(mzp)*hrw(ic)+dhrwp(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsp(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdgsp(mzp)*hrs(ic)+dhrsp(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwp(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdgwp(mzp)*hrw(ic)+dhrwp(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k)) - tzk(mzp)*dtp(ic)
           ablkhh = ablkhh - tz(mzp)*((dzvdsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdsh(mzp)*hrs(ic)+dhrsh(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdwh(mzp)*hrw(ic)+dhrwh(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdgsh(mzp)*hrs(ic)+dhrsh(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdgwh(mzp)*hrw(ic)+dhrwh(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k)) - tzk(mzp)*dth(ic)
           rhsh = rhsh - tz(mzp)*  &
                (((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
                (p(iczp)-p(ic))+  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))  &
                - tzk(mzp)*(tc(iczp)-tc(ic))
        END IF
        IF(ci(1,ma) > 0) THEN
           ablkhp = ablkhp - tz(mzm)*  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm)+  &
                ((dzvdsp(mzm)*hrs(ic)+dhrsp(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsp(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwp(mzm)*hrw(ic)+dhrwp(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwp(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsp(mzm)*hrs(ic)+dhrsp(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsp(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwp(mzm)*hrw(ic)+dhrwp(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwp(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))  &
                - tzk(mzm)*dtp(ic)
           ablkhh = ablkhh - tz(mzm)*(((dzvdsh(mzm)*hrs(ic)+dhrsh(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwh(mzm)*hrw(ic)+dhrwh(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsh(mzm)*hrs(ic)+dhrsh(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwh(mzm)*hrw(ic)+dhrwh(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))  &
                - tzk(mzm)*dth(ic)
           rhsh = rhsh - tz(mzm)*  &
                (((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
                (p(iczm)-p(ic))-  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))  &
                - tzk(mzm)*(tc(iczm)-tc(ic))
        END IF
        ! ... Get seepage boundary face flow rate from residual and calculate advective heat
        ! ...     loss term for heat equation. It is lagged one NR iteration.
        CALL gov_seep(gveqnm,gveqne, i, j, k)
        ! ... Add in the capacitance terms. No source for seeping cell
        gveqnm = gveqnm - (xm(ic)-xmoldt(ic))/delt
        qseepbc = -gveqnm     ! ... Residual net inflow is seepage outflow
        ! ... Calculate the heat sink term for a seeping node.
        ! ... Only the case of single-phase seepage is allowed.
        ! ...      Compressed water or air-water
        ! ... No other heat sources are allowed in the seepage cell
        qhseep = qseepbc*h(ic)
        dqhhseep = qseepbc     ! ... Jacobian terms for heat eq.
        dqhpseep = 0._kdp        ! ... Constant pressure at seeping cell
        vol = dx(i)*dy(j)*dz(k)
        IF(irad) vol = (rsq(i+1)-rsq(i))*pi*dz(k)
        ablkhp = ablkhp - c(ic)*vol/delt + dqhpseep
        ablkhh = ablkhh - d(ic)*vol/delt + dqhhseep
        ! ... Add in the capacitance and b.c. source terms
        rhsh = rhsh + (en(ic)-enoldt(ic))/delt - qhseep
        !    ******* Off Diagonal Block Terms  ***************
        ! ... derivatives of residual flow and transport equations w.r.t. P & H at l+1
        IF(ci(4,ma) > 0) THEN
           ablkpp4 = 0._kdp
           ablkph4 = 0._kdp
           ablkhp4 = tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                ((dxvdsp(mxp)*hrs(icxp)+dhrsp(icxp)*xvds(mxp))*usx(mxp)+  &
                dxvdsp(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                (dxvdwp(mxp)*hrw(icxp)+dhrwp(icxp)*xvdw(mxp))*uwx(mxp)+  &
                dxvdwp(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic)))  &
                + txk(mxp)*dtp(icxp)
           ablkhh4 = tx(mxp)*  &
                ((dxvdsh(mxp)*hrs(icxp)+dhrsh(icxp)*xvds(mxp))*usx(mxp)+  &
                dxvdsh(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                (dxvdwh(mxp)*hrw(icxp)+dhrwh(icxp)*xvdw(mxp))*uwx(mxp)+  &
                dxvdwh(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic))  &
                + txk(mxp)*dth(icxp)
        END IF
        IF(ci(3,ma) > 0) THEN
           ablkpp3 = 0._kdp
           ablkph3 = 0._kdp
           ablkhp3 = tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                (dxvdsp(mxm)*hrs(ic)*usx(mxm)+  &
                (dxvdsp(mxm)*hrs(icxm)+dhrsp(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                dxvdwp(mxm)*hrw(ic)*uwx(mxm)+  &
                (dxvdwp(mxm)*hrw(icxm)+dhrwp(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                (p(icxm)-p(ic))) + txk(mxm)*dtp(icxm)
           ablkhh3 = tx(mxm)*(dxvdsh(mxm)*hrs(ic)*usx(mxm)+  &
                (dxvdsh(mxm)*hrs(icxm)+dhrsh(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                dxvdwh(mxm)*hrw(ic)*uwx(mxm)+  &
                (dxvdwh(mxm)*hrw(icxm)+dhrwh(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                (p(icxm)-p(ic)) + txk(mxm)*dth(icxm)
        END IF
        IF(ci(5,ma) > 0) THEN
           ablkpp5 = 0._kdp
           ablkph5 = 0._kdp
           ablkhp5 = ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                ((dyvdsp(myp)*hrs(icyp)+dhrsp(icyp)*yvds(myp))*usy(myp)+  &
                dyvdsp(myp)*hrs(ic)*(1._kdp-usy(myp))+  &
                (dyvdwp(myp)*hrw(icyp)+dhrwp(icyp)*yvdw(myp))*uwy(myp)+  &
                dyvdwp(myp)*hrw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic)))  &
                + tyk(myp)*dtp(icyp)
           ablkhh5 = ty(myp)*  &
                ((dyvdsh(myp)*hrs(icyp)+dhrsh(icyp)*yvds(myp))*usy(myp)+  &
                dyvdsh(myp)*hrs(ic)*(1._kdp-usy(myp))+  &
                (dyvdwh(myp)*hrw(icyp)+dhrwh(icyp)*yvdw(myp))*uwy(myp)+  &
                dyvdwh(myp)*hrw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic))  &
                + tyk(myp)*dth(icyp)
        END IF
        IF(ci(2,ma) > 0) THEN
           ablkpp2 = 0._kdp
           ablkph2 = 0._kdp
           ablkhp2 = ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                (dyvdsp(mym)*hrs(ic)*usy(mym)+  &
                (dyvdsp(mym)*hrs(icym)+dhrsp(icym)*yvds(mym))*(1._kdp-usy(mym))+  &
                dyvdwp(mym)*hrw(ic)*uwy(mym)+  &
                (dyvdwp(mym)*hrw(icym)+dhrwp(icym)*yvdw(mym))*(1._kdp-uwy(mym)))*  &
                (p(icym)-p(ic))) + tyk(mym)*dtp(icym)
           ablkhh2 = ty(mym)*(dyvdsh(mym)*hrs(ic)*usy(mym)+  &
                (dyvdsh(mym)*hrs(icym)+dhrsh(icym)*yvds(mym))*(1._kdp-usy(mym))+  &
                dyvdwh(mym)*hrw(ic)*uwy(mym)+  &
                (dyvdwh(mym)*hrw(icym)+dhrwh(icym)*yvdw(mym))*(1._kdp-uwy(mym)))*  &
                (p(icym)-p(ic)) + tyk(mym)*dth(icym)
        END IF
        IF(ci(6,ma) > 0) THEN
           ablkpp6 = 0._kdp
           ablkph6 = 0._kdp
           ablkhp6 = tz(mzp)*  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp)+  &
                ((dzvdsp(mzp)*hrs(iczp)+dhrsp(iczp)*zvds(mzp))*usz(mzp)+  &
                dzvdsp(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdwp(mzp)*hrw(iczp)+dhrwp(iczp)*zvdw(mzp))*uwz(mzp)+  &
                dzvdwp(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                ((dzvdgsp(mzp)*hrs(iczp)+dhrsp(iczp)*zvdgs(mzp))*usz(mzp)+  &
                dzvdgsp(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdgwp(mzp)*hrw(iczp)+dhrwp(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                dzvdgwp(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*zzz(k))  &
                + tzk(mzp)*dtp(iczp)
           ablkhh6 = tz(mzp)*  &
                (((dzvdsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvds(mzp))*usz(mzp)+  &
                dzvdsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdw(mzp))*uwz(mzp)+  &
                dzvdwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                ((dzvdgsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvdgs(mzp))*usz(mzp)+  &
                dzvdgsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdgwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                dzvdgwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*zzz(k))  &
                + tzk(mzp)*dth(iczp)
        END IF
        IF(ci(1,ma) > 0) THEN
           ablkpp1 = 0._kdp
           ablkph1 = 0._kdp
           ablkhp1 = tz(mzm)*  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm)+  &
                (dzvdsp(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdsp(mzm)*hrs(iczm)+dhrsp(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                dzvdwp(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdwp(mzm)*hrw(iczm)+dhrwp(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                (p(iczm)-p(ic))- (dzvdgsp(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdgsp(mzm)*hrs(iczm)+dhrsp(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                dzvdgwp(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdgwp(mzm)*hrw(iczm)+dhrwp(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                zzz(k-1)) + tzk(mzm)*dtp(iczm)
           ablkhh1 = tz(mzm)*((dzvdsh(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                dzvdwh(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                (p(iczm)-p(ic))- (dzvdgsh(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdgsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                dzvdgwh(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdgwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                zzz(k-1)) + tzk(mzm)*dth(iczm)
        END IF
     ELSE
        ! ... Normal assembly of active cell
        !    ******* Main Diagonal Block Terms and Right Hand Side *********
        ! ... main diagonal stored first in each row
        IF(ci(4,ma) > 0) THEN
           ! ... Flow
           ablkpp = ablkpp - tx(mxp)*  &
                ((rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                (dxvdsp(mxp)*rs(icxp)*usx(mxp)+  &
                (dxvdsp(mxp)*rs(ic)+drsp(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwp(mxp)*rw(icxp)*uwx(mxp)+  &
                (dxvdwp(mxp)*rw(ic)+drwp(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp)))
           ablkph = ablkph - tx(mxp)*(dxvdsh(mxp)*rs(icxp)*usx(mxp)+  &
                (dxvdsh(mxp)*rs(ic)+drsh(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwh(mxp)*rw(icxp)*uwx(mxp)+  &
                (dxvdwh(mxp)*rw(ic)+drwh(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp))
           ! ... Heat
           ablkhp = ablkhp - tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                (dxvdsp(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsp(mxp)*hrs(ic)+dhrsp(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwp(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwp(mxp)*hrw(ic)+dhrwp(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp))) - txk(mxp)*dtp(ic)
           ablkhh = ablkhh - tx(mxp)*(dxvdsh(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsh(mxp)*hrs(ic)+dhrsh(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwh(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwh(mxp)*hrw(ic)+dhrwh(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp)) - txk(mxp)*dth(ic)
           ! ... Flow
           rhsp = rhsp - tx(mxp)*  &
                ((rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*(p(icxp)-p(ic))
           ! ... Heat
           rhsh = rhsh - tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*  &
                (p(icxp)-p(ic)) - txk(mxp)*(tc(icxp)-tc(ic))
        END IF
        IF(ci(3,ma) > 0) THEN
           ablkpp = ablkpp - tx(mxm)*  &
                ((rs(ic)*usx(mxm)+rs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (rw(ic)*uwx(mxm)+rw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                ((dxvdsp(mxm)*rs(ic)+drsp(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsp(mxm)*rs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwp(mxm)*rw(ic)+drwp(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwp(mxm)*rw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm)))
           ablkph = ablkph - tx(mxm)*((dxvdsh(mxm)*rs(ic)+drsh(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsh(mxm)*rs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwh(mxm)*rw(ic)+drwh(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwh(mxm)*rw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm))
           ablkhp = ablkhp - tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                ((dxvdsp(mxm)*hrs(ic)+dhrsp(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsp(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwp(mxm)*hrw(ic)+dhrwp(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwp(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm)))  &
                - txk(mxm)*dtp(ic)
           ablkhh = ablkhh - tx(mxm)*((dxvdsh(mxm)*hrs(ic)+dhrsh(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsh(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwh(mxm)*hrw(ic)+dhrwh(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwh(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm))  &
                - txk(mxm)*dth(ic)
           rhsp = rhsp - tx(mxm)*  &
                ((rs(ic)*usx(mxm)+rs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (rw(ic)*uwx(mxm)+rw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*(p(icxm)-p(ic))
           rhsh = rhsh - tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*  &
                (p(icxm)-p(ic)) - txk(mxm)*(tc(icxm)-tc(ic))
        END IF
        IF(ci(5,ma) > 0) THEN
           ablkpp = ablkpp - ty(myp)*  &
                ((rs(icyp)*usy(myp)+rs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (rw(icyp)*uwy(myp)+rw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                (dyvdsp(myp)*rs(icyp)*usy(myp)+  &
                (dyvdsp(myp)*rs(ic)+drsp(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwp(myp)*rw(icyp)*uwy(myp)+  &
                (dyvdwp(myp)*rw(ic)+drwp(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp)))
           ablkph = ablkph - ty(myp)*(dyvdsh(myp)*rs(icyp)*usy(myp)+  &
                (dyvdsh(myp)*rs(ic)+drsh(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwh(myp)*rw(icyp)*uwy(myp)+  &
                (dyvdwh(myp)*rw(ic)+drwh(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp))
           ablkhp = ablkhp - ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                (dyvdsp(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsp(myp)*hrs(ic)+dhrsp(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwp(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwp(myp)*hrw(ic)+dhrwp(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp))) - tyk(myp)*dtp(ic)
           ablkhh = ablkhh - ty(myp)*(dyvdsh(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsh(myp)*hrs(ic)+dhrsh(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwh(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwh(myp)*hrw(ic)+dhrwh(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp)) - tyk(myp)*dth(ic)
           rhsp = rhsp - ty(myp)*  &
                ((rs(icyp)*usy(myp)+rs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (rw(icyp)*uwy(myp)+rw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*(p(icyp)-p(ic))
           rhsh = rhsh - ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp))*  &
                (p(icyp)-p(ic)) - tyk(myp)*(tc(icyp)-tc(ic))
        END IF
        IF(ci(2,ma) > 0) THEN
           ablkpp = ablkpp - ty(mym)*  &
                ((rs(ic)*usy(mym)+rs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (rw(ic)*uwy(mym)+rw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                ((dyvdsp(mym)*rs(ic)+drsp(ic)*yvds(mym))*usy(mym)+  &
                dyvdsp(mym)*rs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwp(mym)*rw(ic)+drwp(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwp(mym)*rw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym)))
           ablkph = ablkph - ty(mym)*((dyvdsh(mym)*rs(ic)+drsh(ic)*yvds(mym))*usy(mym)+  &
                dyvdsh(mym)*rs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwh(mym)*rw(ic)+drwh(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwh(mym)*rw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym))
           ablkhp = ablkhp - ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                ((dyvdsp(mym)*hrs(ic)+dhrsp(ic)*yvds(mym))*usy(mym)+  &
                dyvdsp(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwp(mym)*hrw(ic)+dhrwp(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwp(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym)))  &
                - tyk(mym)*dtp(ic)
           ablkhh = ablkhh - ty(mym)*((dyvdsh(mym)*hrs(ic)+dhrsh(ic)*yvds(mym))*usy(mym)+  &
                dyvdsh(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwh(mym)*hrw(ic)+dhrwh(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwh(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym))  &
                - tyk(mym)*dth(ic)
           rhsp = rhsp - ty(mym)*  &
                ((rs(ic)*usy(mym)+rs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (rw(ic)*uwy(mym)+rw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*(p(icym)-p(ic))
           rhsh = rhsh -ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym))*  &
                (p(icym)-p(ic)) - tyk(mym)*(tc(icym)-tc(ic))
        END IF
        IF(ci(6,ma) > 0) THEN
           ablkpp = ablkpp - tz(mzp)*  &
                ((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp)+  &
                (dzvdsp(mzp)*rs(iczp)*usz(mzp)+  &
                (dzvdsp(mzp)*rs(ic)+drsp(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwp(mzp)*rw(iczp)*uwz(mzp)+  &
                (dzvdwp(mzp)*rw(ic)+drwp(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsp(mzp)*rs(iczp)*usz(mzp)+  &
                (dzvdgsp(mzp)*rs(ic)+drsp(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwp(mzp)*rw(iczp)*uwz(mzp)+  &
                (dzvdgwp(mzp)*rw(ic)+drwp(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k))
           ablkph = ablkph - tz(mzp)*((dzvdsh(mzp)*rs(iczp)*usz(mzp)+  &
                (dzvdsh(mzp)*rs(ic)+drsh(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwh(mzp)*rw(iczp)*uwz(mzp)+  &
                (dzvdwh(mzp)*rw(ic)+drwh(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsh(mzp)*rs(iczp)*usz(mzp)+  &
                (dzvdgsh(mzp)*rs(ic)+drsh(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwh(mzp)*rw(iczp)*uwz(mzp)+  &
                (dzvdgwh(mzp)*rw(ic)+drwh(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k))
           ablkhp = ablkhp - tz(mzp)*  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp)+  &
                (dzvdsp(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdsp(mzp)*hrs(ic)+dhrsp(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwp(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdwp(mzp)*hrw(ic)+dhrwp(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsp(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdgsp(mzp)*hrs(ic)+dhrsp(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwp(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdgwp(mzp)*hrw(ic)+dhrwp(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k)) - tzk(mzp)*dtp(ic)
           ablkhh = ablkhh - tz(mzp)*((dzvdsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdsh(mzp)*hrs(ic)+dhrsh(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdwh(mzp)*hrw(ic)+dhrwh(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdgsh(mzp)*hrs(ic)+dhrsh(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdgwh(mzp)*hrw(ic)+dhrwh(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k)) - tzk(mzp)*dth(ic)
           rhsp = rhsp - tz(mzp)*  &
                (((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
                (p(iczp)-p(ic))+  &
                ((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
                (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))
           rhsh = rhsh - tz(mzp)*  &
                (((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
                (p(iczp)-p(ic))+  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))  &
                - tzk(mzp)*(tc(iczp)-tc(ic))
        END IF
        IF(ci(1,ma) > 0) THEN
           ablkpp = ablkpp - tz(mzm)*  &
                ((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm)+  &
                ((dzvdsp(mzm)*rs(ic)+drsp(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsp(mzm)*rs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwp(mzm)*rw(ic)+drwp(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwp(mzm)*rw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsp(mzm)*rs(ic)+drsp(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsp(mzm)*rs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwp(mzm)*rw(ic)+drwp(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwp(mzm)*rw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))
           ablkph = ablkph - tz(mzm)*(((dzvdsh(mzm)*rs(ic)+drsh(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsh(mzm)*rs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwh(mzm)*rw(ic)+drwh(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwh(mzm)*rw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsh(mzm)*rs(ic)+drsh(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsh(mzm)*rs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwh(mzm)*rw(ic)+drwh(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwh(mzm)*rw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))
           ablkhp = ablkhp - tz(mzm)*  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm)+  &
                ((dzvdsp(mzm)*hrs(ic)+dhrsp(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsp(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwp(mzm)*hrw(ic)+dhrwp(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwp(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsp(mzm)*hrs(ic)+dhrsp(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsp(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwp(mzm)*hrw(ic)+dhrwp(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwp(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))  &
                - tzk(mzm)*dtp(ic)
           ablkhh = ablkhh - tz(mzm)*(((dzvdsh(mzm)*hrs(ic)+dhrsh(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwh(mzm)*hrw(ic)+dhrwh(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsh(mzm)*hrs(ic)+dhrsh(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwh(mzm)*hrw(ic)+dhrwh(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))  &
                - tzk(mzm)*dth(ic)
           rhsp = rhsp - tz(mzm)*  &
                (((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
                (p(iczm)-p(ic))-  &
                ((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
                (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))
           rhsh = rhsh - tz(mzm)*  &
                (((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
                (p(iczm)-p(ic))-  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))  &
                - tzk(mzm)*(tc(iczm)-tc(ic))
        END IF
        vol = dx(i)*dy(j)*dz(k)
        IF(irad) vol = (rsq(i+1)-rsq(i))*pi*dz(k)
        ablkpp = ablkpp - f(ic)*vol/delt
        ablkph = ablkph - g(ic)*vol/delt
        ablkhp = ablkhp - c(ic)*vol/delt + dqhp(ic)
        ablkhh = ablkhh - d(ic)*vol/delt + dqhh(ic)
        ! ... Add in the capacitance and source terms
        rhsp = rhsp + (xm(ic)-xmoldt(ic))/delt - q(ic)
        rhsh = rhsh + (en(ic)-enoldt(ic))/delt - qh(ic)
        !    ******* Off Diagonal Block Terms  ***************
        ! ... derivatives of residual flow and transport equations w.r.t. P & H at l+1
        IF(ci(4,ma) > 0) THEN
           ablkpp4 = tx(mxp)*  &
                ((rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                ((dxvdsp(mxp)*rs(icxp)+drsp(icxp)*xvds(mxp))*usx(mxp)+  &
                dxvdsp(mxp)*rs(ic)*(1._kdp-usx(mxp))+  &
                (dxvdwp(mxp)*rw(icxp)+drwp(icxp)*xvdw(mxp))*uwx(mxp)+  &
                dxvdwp(mxp)*rw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic)))
           ablkph4 = tx(mxp)*  &
                ((dxvdsh(mxp)*rs(icxp)+drsh(icxp)*xvds(mxp))*usx(mxp)+  &
                dxvdsh(mxp)*rs(ic)*(1._kdp-usx(mxp))+  &
                (dxvdwh(mxp)*rw(icxp)+drwh(icxp)*xvdw(mxp))*uwx(mxp)+  &
                dxvdwh(mxp)*rw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic))
           ablkhp4 = tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                ((dxvdsp(mxp)*hrs(icxp)+dhrsp(icxp)*xvds(mxp))*usx(mxp)+  &
                dxvdsp(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                (dxvdwp(mxp)*hrw(icxp)+dhrwp(icxp)*xvdw(mxp))*uwx(mxp)+  &
                dxvdwp(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic)))  &
                + txk(mxp)*dtp(icxp)
           ablkhh4 = tx(mxp)*  &
                ((dxvdsh(mxp)*hrs(icxp)+dhrsh(icxp)*xvds(mxp))*usx(mxp)+  &
                dxvdsh(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                (dxvdwh(mxp)*hrw(icxp)+dhrwh(icxp)*xvdw(mxp))*uwx(mxp)+  &
                dxvdwh(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic))  &
                + txk(mxp)*dth(icxp)
        END IF
        IF(ci(3,ma) > 0) THEN
           ablkpp3 = tx(mxm)*  &
                ((rs(ic)*usx(mxm)+rs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (rw(ic)*uwx(mxm)+rw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                (dxvdsp(mxm)*rs(ic)*usx(mxm)+  &
                (dxvdsp(mxm)*rs(icxm)+drsp(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                dxvdwp(mxm)*rw(ic)*uwx(mxm)+  &
                (dxvdwp(mxm)*rw(icxm)+drwp(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                (p(icxm)-p(ic)))
           ablkph3 = tx(mxm)*(dxvdsh(mxm)*rs(ic)*usx(mxm)+  &
                (dxvdsh(mxm)*rs(icxm)+drsh(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                dxvdwh(mxm)*rw(ic)*uwx(mxm)+  &
                (dxvdwh(mxm)*rw(icxm)+drwh(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                (p(icxm)-p(ic))
           ablkhp3 = tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                (dxvdsp(mxm)*hrs(ic)*usx(mxm)+  &
                (dxvdsp(mxm)*hrs(icxm)+dhrsp(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                dxvdwp(mxm)*hrw(ic)*uwx(mxm)+  &
                (dxvdwp(mxm)*hrw(icxm)+dhrwp(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                (p(icxm)-p(ic))) + txk(mxm)*dtp(icxm)
           ablkhh3 = tx(mxm)*(dxvdsh(mxm)*hrs(ic)*usx(mxm)+  &
                (dxvdsh(mxm)*hrs(icxm)+dhrsh(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                dxvdwh(mxm)*hrw(ic)*uwx(mxm)+  &
                (dxvdwh(mxm)*hrw(icxm)+dhrwh(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                (p(icxm)-p(ic)) + txk(mxm)*dth(icxm)
        END IF
        IF(ci(5,ma) > 0) THEN
           ablkpp5 = ty(myp)*  &
                ((rs(icyp)*usy(myp)+rs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (rw(icyp)*uwy(myp)+rw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                ((dyvdsp(myp)*rs(icyp)+drsp(icyp)*yvds(myp))*usy(myp)+  &
                dyvdsp(myp)*rs(ic)*(1._kdp-usy(myp))+  &
                (dyvdwp(myp)*rw(icyp)+drwp(icyp)*yvdw(myp))*uwy(myp)+  &
                dyvdwp(myp)*rw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic)))
           ablkph5 = ty(myp)*  &
                ((dyvdsh(myp)*rs(icyp)+drsh(icyp)*yvds(myp))*usy(myp)+  &
                dyvdsh(myp)*rs(ic)*(1._kdp-usy(myp))+  &
                (dyvdwh(myp)*rw(icyp)+drwh(icyp)*yvdw(myp))*uwy(myp)+  &
                dyvdwh(myp)*rw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic))
           ablkhp5 = ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                ((dyvdsp(myp)*hrs(icyp)+dhrsp(icyp)*yvds(myp))*usy(myp)+  &
                dyvdsp(myp)*hrs(ic)*(1._kdp-usy(myp))+  &
                (dyvdwp(myp)*hrw(icyp)+dhrwp(icyp)*yvdw(myp))*uwy(myp)+  &
                dyvdwp(myp)*hrw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic)))  &
                + tyk(myp)*dtp(icyp)
           ablkhh5 = ty(myp)*  &
                ((dyvdsh(myp)*hrs(icyp)+dhrsh(icyp)*yvds(myp))*usy(myp)+  &
                dyvdsh(myp)*hrs(ic)*(1._kdp-usy(myp))+  &
                (dyvdwh(myp)*hrw(icyp)+dhrwh(icyp)*yvdw(myp))*uwy(myp)+  &
                dyvdwh(myp)*hrw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic))  &
                + tyk(myp)*dth(icyp)
        END IF
        IF(ci(2,ma) > 0) THEN
           ablkpp2 = ty(mym)*  &
                ((rs(ic)*usy(mym)+rs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (rw(ic)*uwy(mym)+rw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                (dyvdsp(mym)*rs(ic)*usy(mym)+  &
                (dyvdsp(mym)*rs(icym)+drsp(icym)*yvds(mym))*(1._kdp-usy(mym))+  &
                dyvdwp(mym)*rw(ic)*uwy(mym)+  &
                (dyvdwp(mym)*rw(icym)+drwp(icym)*yvdw(mym))*(1._kdp-uwy(mym)))*  &
                (p(icym)-p(ic)))
           ablkph2 = ty(mym)*(dyvdsh(mym)*rs(ic)*usy(mym)+  &
                (dyvdsh(mym)*rs(icym)+drsh(icym)*yvds(mym))*(1._kdp-usy(mym))+  &
                dyvdwh(mym)*rw(ic)*uwy(mym)+  &
                (dyvdwh(mym)*rw(icym)+drwh(icym)*yvdw(mym))*(1._kdp-uwy(mym)))*  &
                (p(icym)-p(ic))
           ablkhp2 = ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                (dyvdsp(mym)*hrs(ic)*usy(mym)+  &
                (dyvdsp(mym)*hrs(icym)+dhrsp(icym)*yvds(mym))*(1._kdp-usy(mym))+  &
                dyvdwp(mym)*hrw(ic)*uwy(mym)+  &
                (dyvdwp(mym)*hrw(icym)+dhrwp(icym)*yvdw(mym))*(1._kdp-uwy(mym)))*  &
                (p(icym)-p(ic))) + tyk(mym)*dtp(icym)
           ablkhh2 = ty(mym)*(dyvdsh(mym)*hrs(ic)*usy(mym)+  &
                (dyvdsh(mym)*hrs(icym)+dhrsh(icym)*yvds(mym))*(1._kdp-usy(mym))+  &
                dyvdwh(mym)*hrw(ic)*uwy(mym)+  &
                (dyvdwh(mym)*hrw(icym)+dhrwh(icym)*yvdw(mym))*(1._kdp-uwy(mym)))*  &
                (p(icym)-p(ic)) + tyk(mym)*dth(icym)
        END IF
        IF(ci(6,ma) > 0) THEN
           ablkpp6 = tz(mzp)*  &
                ((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp)+  &
                ((dzvdsp(mzp)*rs(iczp)+drsp(iczp)*zvds(mzp))*usz(mzp)+  &
                dzvdsp(mzp)*rs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdwp(mzp)*rw(iczp)+drwp(iczp)*zvdw(mzp))*uwz(mzp)+  &
                dzvdwp(mzp)*rw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                ((dzvdgsp(mzp)*rs(iczp)+drsp(iczp)*zvdgs(mzp))*usz(mzp)+  &
                dzvdgsp(mzp)*rs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdgwp(mzp)*rw(iczp)+drwp(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                dzvdgwp(mzp)*rw(ic)*(1._kdp-uwz(mzp)))*zzz(k))
           ablkph6 = tz(mzp)*  &
                (((dzvdsh(mzp)*rs(iczp)+drsh(iczp)*zvds(mzp))*usz(mzp)+  &
                dzvdsh(mzp)*rs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdwh(mzp)*rw(iczp)+drwh(iczp)*zvdw(mzp))*uwz(mzp)+  &
                dzvdwh(mzp)*rw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                ((dzvdgsh(mzp)*rs(iczp)+drsh(iczp)*zvdgs(mzp))*usz(mzp)+  &
                dzvdgsh(mzp)*rs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdgwh(mzp)*rw(iczp)+drwh(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                dzvdgwh(mzp)*rw(ic)*(1._kdp-uwz(mzp)))*zzz(k))
           ablkhp6 = tz(mzp)*  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp)+  &
                ((dzvdsp(mzp)*hrs(iczp)+dhrsp(iczp)*zvds(mzp))*usz(mzp)+  &
                dzvdsp(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdwp(mzp)*hrw(iczp)+dhrwp(iczp)*zvdw(mzp))*uwz(mzp)+  &
                dzvdwp(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                ((dzvdgsp(mzp)*hrs(iczp)+dhrsp(iczp)*zvdgs(mzp))*usz(mzp)+  &
                dzvdgsp(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdgwp(mzp)*hrw(iczp)+dhrwp(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                dzvdgwp(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*zzz(k))  &
                + tzk(mzp)*dtp(iczp)
           ablkhh6 = tz(mzp)*  &
                (((dzvdsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvds(mzp))*usz(mzp)+  &
                dzvdsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdw(mzp))*uwz(mzp)+  &
                dzvdwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                ((dzvdgsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvdgs(mzp))*usz(mzp)+  &
                dzvdgsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                (dzvdgwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                dzvdgwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*zzz(k))  &
                + tzk(mzp)*dth(iczp)
        END IF
        IF(ci(1,ma) > 0) THEN
           ablkpp1 = tz(mzm)*  &
                ((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm)+  &
                (dzvdsp(mzm)*rs(ic)*usz(mzm)+  &
                (dzvdsp(mzm)*rs(iczm)+drsp(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                dzvdwp(mzm)*rw(ic)*uwz(mzm)+  &
                (dzvdwp(mzm)*rw(iczm)+drwp(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                (p(iczm)-p(ic))- (dzvdgsp(mzm)*rs(ic)*usz(mzm)+  &
                (dzvdgsp(mzm)*rs(iczm)+drsp(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                dzvdgwp(mzm)*rw(ic)*uwz(mzm)+  &
                (dzvdgwp(mzm)*rw(iczm)+drwp(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                zzz(k-1))
           ablkph1 = tz(mzm)*((dzvdsh(mzm)*rs(ic)*usz(mzm)+  &
                (dzvdsh(mzm)*rs(iczm)+drsh(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                dzvdwh(mzm)*rw(ic)*uwz(mzm)+  &
                (dzvdwh(mzm)*rw(iczm)+drwh(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                (p(iczm)-p(ic))- (dzvdgsh(mzm)*rs(ic)*usz(mzm)+  &
                (dzvdgsh(mzm)*rs(iczm)+drsh(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                dzvdgwh(mzm)*rw(ic)*uwz(mzm)+  &
                (dzvdgwh(mzm)*rw(iczm)+drwh(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                zzz(k-1))
           ablkhp1 = tz(mzm)*  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm)+  &
                (dzvdsp(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdsp(mzm)*hrs(iczm)+dhrsp(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                dzvdwp(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdwp(mzm)*hrw(iczm)+dhrwp(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                (p(iczm)-p(ic))- (dzvdgsp(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdgsp(mzm)*hrs(iczm)+dhrsp(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                dzvdgwp(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdgwp(mzm)*hrw(iczm)+dhrwp(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                zzz(k-1)) + tzk(mzm)*dtp(iczm)
           ablkhh1 = tz(mzm)*((dzvdsh(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                dzvdwh(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                (p(iczm)-p(ic))- (dzvdgsh(mzm)*hrs(ic)*usz(mzm)+  &
                (dzvdgsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                dzvdgwh(mzm)*hrw(ic)*uwz(mzm)+  &
                (dzvdgwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                zzz(k-1)) + tzk(mzm)*dth(iczm)
        END IF
     END IF
     ! ... Assembly of active cell; pressure
     ! ...      Load rhs array
     ! ... Make diagonal elements positive by multiplying equation by -1
     rhs(2*ma-1) = -rhsp
     ! ... Load values into av array
     ! ... Diagonal array elements
     ! ... main diagonal stored first in each row
     nnn=nnn+1
     av(nnn) = -ablkpp
     nnn=nnn+1
     av(nnn) = -ablkph
     ! ...  Off-Diagonal array elements
     IF(ci(1,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkpp1
        nnn=nnn+1
        av(nnn) = -ablkph1
     END IF
     IF(ci(2,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkpp2
        nnn=nnn+1
        av(nnn) = -ablkph2
     END IF
     IF(ci(3,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkpp3
        nnn=nnn+1
        av(nnn) = -ablkph3
     END IF
     IF(ci(4,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkpp4
        nnn=nnn+1
        av(nnn) = -ablkph4
     END IF
     IF(ci(5,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkpp5
        nnn=nnn+1
        av(nnn) = -ablkph5
     END IF
     IF(ci(6,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkpp6
        nnn=nnn+1
        av(nnn) = -ablkph6
     END IF
     ! ...  Assembly of active cell; enthalpy
     rhs(2*ma) = -rhsh
     ! ... Diagonal array elements
     nnn=nnn+1
     av(nnn) = -ablkhh      ! ... diagonal is loaded first
     nnn=nnn+1
     av(nnn) = -ablkhp
     ! ...  Off-Diagonal array elements
     IF(ci(1,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkhp1
        nnn=nnn+1
        av(nnn) = -ablkhh1
     END IF
     IF(ci(2,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkhp2
        nnn=nnn+1
        av(nnn) = -ablkhh2
     END IF
     IF(ci(3,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkhp3
        nnn=nnn+1
        av(nnn) = -ablkhh3
     END IF
     IF(ci(4,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkhp4
        nnn=nnn+1
        av(nnn) = -ablkhh4
     END IF
     IF(ci(5,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkhp5
        nnn=nnn+1
        av(nnn) = -ablkhh5
     END IF
     IF(ci(6,ma) > 0) THEN
        nnn=nnn+1
        av(nnn) = -ablkhp6
        nnn=nnn+1
        av(nnn) = -ablkhh6
     END IF
  END DO
END SUBROUTINE assemble_nre_csr
