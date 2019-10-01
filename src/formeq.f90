SUBROUTINE formeq(j)
  !     Purpose: To assemble the set of Newton-Raphson matrix equations.
  !              Based on the fully implicit finite difference
  !              equations for flow and energy transport and the four
  !              partial derivative equations (partials of the residual flow 
  !              and transport equations with respect to pressure and enthalpy).
  !  ***  use sparskit scaling, do someday ***
  USE machine_constants, ONLY: kdp
  USE math_constants
  USE bc
  USE control, ONLY: ierr
  USE fdeq
  USE mesh
  USE parameters
  USE solver
  USE source
  USE variables
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: j      ! ... slice j of the mesh
  ! 
  INTEGER :: i, ic, icxm, icxp, icym, icyp, iczm, iczp, ik, ikxm, ikxp,  &
       ikzm, ikzp, k, ml, mm, mp
  INTEGER :: mxm, mxp, mym, myp, mzm, mzp, nby, nc, nn
  INTEGER :: ibb, mmr, nnn
  INTEGER :: ii, jj, kl, ku, info
  INTEGER, DIMENSION(1) :: loc_min
  INTEGER :: a_err
  REAL(KIND=kdp) :: a1, a2, a3, a4, r1, r2, vol, ainv
  REAL(KIND=kdp) :: gveqnm, gveqne, qsppbc, qhspp, dqhhspp, dqhpspp,  &
       qseepbc, qhseep, dqhhseep, dqhpseep 
  REAL(KIND=kdp) :: bignum, smalnum, rcmin, rcmax, aamax, rowcond, colcond, cj
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: r
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.3 $//$Date: 2007/08/24 00:01:59 $'
  !     ------------------------------------------------------------------
  !...
  mp = (mbw+1)/2     ! ...    main diagonal (MP) - stored as array column
  ibb = (mbw-1)/2
  ! ... Initialize arrays R and A
  DO  ml = 1,2*nxz
     rlit(ml) = 0._kdp
     DO  mm = 1,mbw
        alit(ml,mm) = 0._kdp
     END DO
  END DO
!!$  !*****special output for solver testing
!!$      open(60,file='vaht.dat')
  nnn = 0          ! ... counter for array element
  ! ... Assemble for each active node in this slice (NBY)
  nby = nb(j)
  DO  ml = 1,nby
     ! ... Calculate the I and K indices of each active node to be assembled
     i = npp(2*ml-1, j)
     k = npp(2*ml, j)
     ! ... Global node indices
     ic = (k-1)*nx + i + (j-1)*nxz
     ! ... Adjacent node indices
     icxp = ic + 1
     icxm = ic - 1
     iczp = ic + nx
     iczm = ic - nx
     icyp = ic + nxz
     icym = ic - nxz
     ! ... Node indices for current slice
     ik = (k-1)*nx + i
     ikxp = ik + 1
     ikxm = ik - 1
     ikzp = ik + nx
     ikzm = ik - nx
     ! ... Cell face indices for conductance factors
     mxm = (k-1)*nxx + i + (j-1)*nxx*nz
     mxp = mxm + 1
     mym = ic
     myp = mym + nxz
     mzm = (k-1)*nx + i + (j-1)*nxzz
     mzp = mzm + nx
     IF((ibc(ik,j)/10 == 3 .OR. ibc(ik,j)/10 == 5) .AND. seeping(ik,j)) THEN
        ! ... Assembly of seeping cell
        !    ******* Main Diagonal Block Terms and Right Hand Side *********
        ! ... Flow N-R equation
        ! ...      Assemble a trivial equation for pressure 
        alit(2*ml-1,mp) = -1._kdp
        alit(2*ml-1,mp+1) = 0._kdp
        rlit(2*ml-1) = 0._kdp
        ! ... Heat
        IF (i < nx) THEN
           alit(2*ml,mp-1) = -tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                (dxvdsp(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsp(mxp)*hrs(ic)+dhrsp(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwp(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwp(mxp)*hrw(ic)+dhrwp(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp))) - txk(mxp)*dtp(ic)
           alit(2*ml,mp) = -tx(mxp)*(dxvdsh(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsh(mxp)*hrs(ic)+dhrsh(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwh(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwh(mxp)*hrw(ic)+dhrwh(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp)) - txk(mxp)*dth(ic)
           rlit(2*ml) = -tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*  &
                (p(icxp)-p(ic)) - txk(mxp)*(tc(icxp)-tc(ic))
        END IF
        IF (i > 1) THEN
           a3 = -tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                ((dxvdsp(mxm)*hrs(ic)+dhrsp(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsp(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwp(mxm)*hrw(ic)+dhrwp(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwp(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm)))  &
                - txk(mxm)*dtp(ic)
           alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
           a4 = -tx(mxm)*((dxvdsh(mxm)*hrs(ic)+dhrsh(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsh(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwh(mxm)*hrw(ic)+dhrwh(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwh(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm))  &
                - txk(mxm)*dth(ic)
           alit(2*ml,mp) = alit(2*ml,mp) + a4
           r2 = -tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*  &
                (p(icxm)-p(ic)) - txk(mxm)*(tc(icxm)-tc(ic))
           rlit(2*ml) = rlit(2*ml) + r2
        END IF
        IF (j < ny) THEN
           alit(2*ml,mp-1) = alit(2*ml,mp-1) - ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                (dyvdsp(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsp(myp)*hrs(ic)+dhrsp(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwp(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwp(myp)*hrw(ic)+dhrwp(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp))) - tyk(myp)*dtp(ic)
           alit(2*ml,mp) = alit(2*ml,mp) - ty(myp)*  &
                (dyvdsh(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsh(myp)*hrs(ic)+dhrsh(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwh(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwh(myp)*hrw(ic)+dhrwh(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp)) - tyk(myp)*dth(ic)
        END IF
        IF (j > 1) THEN
           alit(2*ml,mp-1) = alit(2*ml,mp-1) - ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                ((dyvdsp(mym)*hrs(ic)+dhrsp(ic)*yvds(mym))*usy(mym)+  &
                dyvdsp(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwp(mym)*hrw(ic)+dhrwp(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwp(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym)))  &
                - tyk(mym)*dtp(ic)
           alit(2*ml,mp) = alit(2*ml,mp) - ty(mym)*  &
                ((dyvdsh(mym)*hrs(ic)+dhrsh(ic)*yvds(mym))*usy(mym)+  &
                dyvdsh(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwh(mym)*hrw(ic)+dhrwh(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwh(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym))  &
                - tyk(mym)*dth(ic)
        END IF
        IF (k < nz) THEN
           a3 = -tz(mzp)*  &
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
           alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
           a4 = -tz(mzp)*((dzvdsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdsh(mzp)*hrs(ic)+dhrsh(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdwh(mzp)*hrw(ic)+dhrwh(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdgsh(mzp)*hrs(ic)+dhrsh(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdgwh(mzp)*hrw(ic)+dhrwh(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k)) - tzk(mzp)*dth(ic)
           alit(2*ml,mp) = alit(2*ml,mp) + a4
           r2 = -tz(mzp)*  &
                (((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
                (p(iczp)-p(ic))+  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))  &
                - tzk(mzp)*(tc(iczp)-tc(ic))
           rlit(2*ml) = rlit(2*ml) + r2
        END IF
        IF (k > 1) THEN
           a3 = -tz(mzm)*  &
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
           alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
           a4 = -tz(mzm)*(((dzvdsh(mzm)*hrs(ic)+dhrsh(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwh(mzm)*hrw(ic)+dhrwh(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsh(mzm)*hrs(ic)+dhrsh(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwh(mzm)*hrw(ic)+dhrwh(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))  &
                - tzk(mzm)*dth(ic)
           alit(2*ml,mp) = alit(2*ml,mp) + a4
           r2 = -tz(mzm)*  &
                (((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
                (p(iczm)-p(ic))-  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))  &
                - tzk(mzm)*(tc(iczm)-tc(ic))
           rlit(2*ml) = rlit(2*ml) + r2
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
        IF (irad) vol = (rsq(i+1)-rsq(i))*pi*dz(k)
        a3 = -c(ic)*vol/delt + dqhpseep
        a4 = -d(ic)*vol/delt + dqhhseep
        alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
        alit(2*ml,mp) = alit(2*ml,mp) + a4
!!$        !... special output
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml-1, 2*ml-1+mp-(ibb+1), Alit(2*ml-1, mp)
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml-1, 2*ml-1+mp+1-(ibb+1), Alit(2*ml-1, mp+1)
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml, 2*ml+mp-1-(ibb+1), Alit(2*ml, mp-1)
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml, 2*ml+mp-(ibb+1), Alit(2*ml, mp)
!!$        !...
        ! ... Add in the capacitance and b.c. source terms
        r2 = (en(ic)-enoldt(ic))/delt - qhseep
        rlit(2*ml) = rlit(2*ml) + r2
        !    ******* Off Diagonal Block Terms  ***************
        ! ... derivatives of residual flow and transport equations w.r.t. P & H at l+1
        IF (i < nx) THEN
           nn = np(ikxp,j)
           IF (nn > 0) THEN
              nc = mp + (nn-ml)*2
              alit(2*ml-1,nc) = 0._kdp
              alit(2*ml-1,nc+1) = 0._kdp
              alit(2*ml,nc-1) = tx(mxp)*  &
                   ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                   (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                   ((dxvdsp(mxp)*hrs(icxp)+dhrsp(icxp)*xvds(mxp))*usx(mxp)+  &
                   dxvdsp(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                   (dxvdwp(mxp)*hrw(icxp)+dhrwp(icxp)*xvdw(mxp))*uwx(mxp)+  &
                   dxvdwp(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic)))  &
                   + txk(mxp)*dtp(icxp)
              alit(2*ml,nc) = tx(mxp)*  &
                   ((dxvdsh(mxp)*hrs(icxp)+dhrsh(icxp)*xvds(mxp))*usx(mxp)+  &
                   dxvdsh(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                   (dxvdwh(mxp)*hrw(icxp)+dhrwh(icxp)*xvdw(mxp))*uwx(mxp)+  &
                   dxvdwh(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic))  &
                   + txk(mxp)*dth(icxp)
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
        IF (i > 1) THEN
           nn = np(ikxm,j)
           IF (nn > 0) THEN
              nc = mp + (nn-ml)*2
              alit(2*ml-1,nc) = 0._kdp
              alit(2*ml-1,nc+1) = 0._kdp
              alit(2*ml,nc-1) = tx(mxm)*  &
                   ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                   (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                   (dxvdsp(mxm)*hrs(ic)*usx(mxm)+  &
                   (dxvdsp(mxm)*hrs(icxm)+dhrsp(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                   dxvdwp(mxm)*hrw(ic)*uwx(mxm)+  &
                   (dxvdwp(mxm)*hrw(icxm)+dhrwp(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                   (p(icxm)-p(ic))) + txk(mxm)*dtp(icxm)
              alit(2*ml,nc) = tx(mxm)*(dxvdsh(mxm)*hrs(ic)*usx(mxm)+  &
                   (dxvdsh(mxm)*hrs(icxm)+dhrsh(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                   dxvdwh(mxm)*hrw(ic)*uwx(mxm)+  &
                   (dxvdwh(mxm)*hrw(icxm)+dhrwh(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                   (p(icxm)-p(ic)) + txk(mxm)*dth(icxm)

!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
        IF (k < nz) THEN
           nn = np(ikzp,j)
           nc = mp + (nn-ml)*2
           IF (nn > 0) THEN
              alit(2*ml-1,nc) = 0._kdp
              alit(2*ml-1,nc+1) = 0._kdp
              alit(2*ml,nc-1) = tz(mzp)*  &
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
              alit(2*ml,nc) = tz(mzp)*  &
                   (((dzvdsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvds(mzp))*usz(mzp)+  &
                   dzvdsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                   (dzvdwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdw(mzp))*uwz(mzp)+  &
                   dzvdwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                   ((dzvdgsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvdgs(mzp))*usz(mzp)+  &
                   dzvdgsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                   (dzvdgwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                   dzvdgwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*zzz(k))  &
                   + tzk(mzp)*dth(iczp)
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
        IF (k > 1) THEN
           nn = np(ikzm,j)
           nc = mp + (nn-ml)*2
           IF (nn > 0) THEN
              alit(2*ml-1,nc) = 0._kdp
              alit(2*ml-1,nc+1) = 0._kdp
              alit(2*ml,nc-1) = tz(mzm)*  &
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
              alit(2*ml,nc) = tz(mzm)*((dzvdsh(mzm)*hrs(ic)*usz(mzm)+  &
                   (dzvdsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                   dzvdwh(mzm)*hrw(ic)*uwz(mzm)+  &
                   (dzvdwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                   (p(iczm)-p(ic))- (dzvdgsh(mzm)*hrs(ic)*usz(mzm)+  &
                   (dzvdgsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                   dzvdgwh(mzm)*hrw(ic)*uwz(mzm)+  &
                   (dzvdgwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                   zzz(k-1)) + tzk(mzm)*dth(iczm)
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
     ELSEIF(ibc(ik,j) == 10) THEN
        ! ... Assembly of specified pressure cell with associated enthalpy
        !    ******* Main Diagonal Block Terms and Right Hand Side *********
        ! ... Flow N-R equation
        ! ...      Assemble a trivial equation for pressure 
        alit(2*ml-1,mp) = -1._kdp
        alit(2*ml-1,mp+1) = 0._kdp
        rlit(2*ml-1) = 0._kdp
        ! ... Heat
        IF (i < nx) THEN
           alit(2*ml,mp-1) = -tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                (dxvdsp(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsp(mxp)*hrs(ic)+dhrsp(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwp(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwp(mxp)*hrw(ic)+dhrwp(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp))) - txk(mxp)*dtp(ic)
           alit(2*ml,mp) = -tx(mxp)*(dxvdsh(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsh(mxp)*hrs(ic)+dhrsh(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwh(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwh(mxp)*hrw(ic)+dhrwh(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp)) - txk(mxp)*dth(ic)
           rlit(2*ml) = -tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*  &
                (p(icxp)-p(ic)) - txk(mxp)*(tc(icxp)-tc(ic))
        END IF
        IF (i > 1) THEN
           a3 = -tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                ((dxvdsp(mxm)*hrs(ic)+dhrsp(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsp(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwp(mxm)*hrw(ic)+dhrwp(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwp(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm)))  &
                - txk(mxm)*dtp(ic)
           alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
           a4 = -tx(mxm)*((dxvdsh(mxm)*hrs(ic)+dhrsh(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsh(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwh(mxm)*hrw(ic)+dhrwh(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwh(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm))  &
                - txk(mxm)*dth(ic)
           alit(2*ml,mp) = alit(2*ml,mp) + a4
           r2 = -tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*  &
                (p(icxm)-p(ic)) - txk(mxm)*(tc(icxm)-tc(ic))
           rlit(2*ml) = rlit(2*ml) + r2
        END IF
        IF (j < ny) THEN
           alit(2*ml,mp-1) = alit(2*ml,mp-1) - ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                (dyvdsp(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsp(myp)*hrs(ic)+dhrsp(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwp(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwp(myp)*hrw(ic)+dhrwp(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp))) - tyk(myp)*dtp(ic)
           alit(2*ml,mp) = alit(2*ml,mp) - ty(myp)*  &
                (dyvdsh(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsh(myp)*hrs(ic)+dhrsh(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwh(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwh(myp)*hrw(ic)+dhrwh(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp)) - tyk(myp)*dth(ic)
        END IF
        IF (j > 1) THEN
           alit(2*ml,mp-1) = alit(2*ml,mp-1) - ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                ((dyvdsp(mym)*hrs(ic)+dhrsp(ic)*yvds(mym))*usy(mym)+  &
                dyvdsp(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwp(mym)*hrw(ic)+dhrwp(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwp(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym)))  &
                - tyk(mym)*dtp(ic)
           alit(2*ml,mp) = alit(2*ml,mp) - ty(mym)*  &
                ((dyvdsh(mym)*hrs(ic)+dhrsh(ic)*yvds(mym))*usy(mym)+  &
                dyvdsh(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwh(mym)*hrw(ic)+dhrwh(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwh(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym))  &
                - tyk(mym)*dth(ic)
        END IF
        IF (k < nz) THEN
           a3 = -tz(mzp)*  &
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
           alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
           a4 = -tz(mzp)*((dzvdsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdsh(mzp)*hrs(ic)+dhrsh(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdwh(mzp)*hrw(ic)+dhrwh(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdgsh(mzp)*hrs(ic)+dhrsh(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdgwh(mzp)*hrw(ic)+dhrwh(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k)) - tzk(mzp)*dth(ic)
           alit(2*ml,mp) = alit(2*ml,mp) + a4
           r2 = -tz(mzp)*  &
                (((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
                (p(iczp)-p(ic))+  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))  &
                - tzk(mzp)*(tc(iczp)-tc(ic))
           rlit(2*ml) = rlit(2*ml) + r2
        END IF
        IF (k > 1) THEN
           a3 = -tz(mzm)*  &
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
           alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
           a4 = -tz(mzm)*(((dzvdsh(mzm)*hrs(ic)+dhrsh(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwh(mzm)*hrw(ic)+dhrwh(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsh(mzm)*hrs(ic)+dhrsh(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwh(mzm)*hrw(ic)+dhrwh(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))  &
                - tzk(mzm)*dth(ic)
           alit(2*ml,mp) = alit(2*ml,mp) + a4
           r2 = -tz(mzm)*  &
                (((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
                (p(iczm)-p(ic))-  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))  &
                - tzk(mzm)*(tc(iczm)-tc(ic))
           rlit(2*ml) = rlit(2*ml) + r2
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
           dqhhspp = qsppbc        ! ... Jacobian terms for heat eq.
        END IF
        dqhpspp = 0._kdp     ! ... Constant pressure at cell so dqhp*delp=0
        vol = dx(i)*dy(j)*dz(k)
        IF (irad) vol = (rsq(i+1)-rsq(i))*pi*dz(k)
        a3 = -c(ic)*vol/delt + dqhpspp
        a4 = -d(ic)*vol/delt + dqhhspp
        alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
        alit(2*ml,mp) = alit(2*ml,mp) + a4
!!$        !... special output
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml-1, 2*ml-1+mp-(ibb+1), Alit(2*ml-1, mp)
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml-1, 2*ml-1+mp+1-(ibb+1), Alit(2*ml-1, mp+1)
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml, 2*ml+mp-1-(ibb+1), Alit(2*ml, mp-1)
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml, 2*ml+mp-(ibb+1), Alit(2*ml, mp)
!!$        !...
        ! ... Add in the capacitance and b.c. source terms
        r2 = (en(ic)-enoldt(ic))/delt - qhspp
        rlit(2*ml) = rlit(2*ml) + r2
        !    ******* Off Diagonal Block Terms  ***************
        ! ... derivatives of residual flow and transport equations w.r.t. P & H at l+1
        IF (i < nx) THEN
           nn = np(ikxp,j)
           IF (nn > 0) THEN
              nc = mp + (nn-ml)*2
              alit(2*ml-1,nc) = 0._kdp
              alit(2*ml-1,nc+1) = 0._kdp
              alit(2*ml,nc-1) = tx(mxp)*  &
                   ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                   (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                   ((dxvdsp(mxp)*hrs(icxp)+dhrsp(icxp)*xvds(mxp))*usx(mxp)+  &
                   dxvdsp(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                   (dxvdwp(mxp)*hrw(icxp)+dhrwp(icxp)*xvdw(mxp))*uwx(mxp)+  &
                   dxvdwp(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic)))  &
                   + txk(mxp)*dtp(icxp)
              alit(2*ml,nc) = tx(mxp)*  &
                   ((dxvdsh(mxp)*hrs(icxp)+dhrsh(icxp)*xvds(mxp))*usx(mxp)+  &
                   dxvdsh(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                   (dxvdwh(mxp)*hrw(icxp)+dhrwh(icxp)*xvdw(mxp))*uwx(mxp)+  &
                   dxvdwh(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic))  &
                   + txk(mxp)*dth(icxp)
!!$
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
        IF (i > 1) THEN
           nn = np(ikxm,j)
           IF (nn > 0) THEN
              nc = mp + (nn-ml)*2
              alit(2*ml-1,nc) = 0._kdp
              alit(2*ml-1,nc+1) = 0._kdp
              alit(2*ml,nc-1) = tx(mxm)*  &
                   ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                   (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                   (dxvdsp(mxm)*hrs(ic)*usx(mxm)+  &
                   (dxvdsp(mxm)*hrs(icxm)+dhrsp(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                   dxvdwp(mxm)*hrw(ic)*uwx(mxm)+  &
                   (dxvdwp(mxm)*hrw(icxm)+dhrwp(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                   (p(icxm)-p(ic))) + txk(mxm)*dtp(icxm)
              alit(2*ml,nc) = tx(mxm)*(dxvdsh(mxm)*hrs(ic)*usx(mxm)+  &
                   (dxvdsh(mxm)*hrs(icxm)+dhrsh(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                   dxvdwh(mxm)*hrw(ic)*uwx(mxm)+  &
                   (dxvdwh(mxm)*hrw(icxm)+dhrwh(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                   (p(icxm)-p(ic)) + txk(mxm)*dth(icxm)
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1,2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
        IF (k < nz) THEN
           nn = np(ikzp,j)
           nc = mp + (nn-ml)*2
           IF (nn > 0) THEN
              alit(2*ml-1,nc) = 0._kdp
              alit(2*ml-1,nc+1) = 0._kdp
              alit(2*ml,nc-1) = tz(mzp)*  &
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
              alit(2*ml,nc) = tz(mzp)*  &
                   (((dzvdsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvds(mzp))*usz(mzp)+  &
                   dzvdsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                   (dzvdwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdw(mzp))*uwz(mzp)+  &
                   dzvdwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                   ((dzvdgsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvdgs(mzp))*usz(mzp)+  &
                   dzvdgsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                   (dzvdgwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                   dzvdgwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*zzz(k))  &
                   + tzk(mzp)*dth(iczp)
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
        IF (k > 1) THEN
           nn = np(ikzm,j)
           nc = mp + (nn-ml)*2
           IF (nn > 0) THEN
              alit(2*ml-1,nc) = 0._kdp
              alit(2*ml-1,nc+1) = 0._kdp
              alit(2*ml,nc-1) = tz(mzm)*  &
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
              alit(2*ml,nc) = tz(mzm)*((dzvdsh(mzm)*hrs(ic)*usz(mzm)+  &
                   (dzvdsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                   dzvdwh(mzm)*hrw(ic)*uwz(mzm)+  &
                   (dzvdwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                   (p(iczm)-p(ic))- (dzvdgsh(mzm)*hrs(ic)*usz(mzm)+  &
                   (dzvdgsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                   dzvdgwh(mzm)*hrw(ic)*uwz(mzm)+  &
                   (dzvdgwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                   zzz(k-1)) + tzk(mzm)*dth(iczm)
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
     ELSE
        ! ... Normal assembly of active cell
        !    ******* Main Diagonal Block Terms and Right Hand Side *********
        IF (i < nx) THEN
           ! ... Flow
           alit(2*ml-1,mp) = -tx(mxp)*  &
                ((rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                (dxvdsp(mxp)*rs(icxp)*usx(mxp)+  &
                (dxvdsp(mxp)*rs(ic)+drsp(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwp(mxp)*rw(icxp)*uwx(mxp)+  &
                (dxvdwp(mxp)*rw(ic)+drwp(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp)))
           alit(2*ml-1,mp+1) = -tx(mxp)*(dxvdsh(mxp)*rs(icxp)*usx(mxp)+  &
                (dxvdsh(mxp)*rs(ic)+drsh(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwh(mxp)*rw(icxp)*uwx(mxp)+  &
                (dxvdwh(mxp)*rw(ic)+drwh(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp))
           ! ... Heat
           alit(2*ml,mp-1) = -tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                (dxvdsp(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsp(mxp)*hrs(ic)+dhrsp(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwp(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwp(mxp)*hrw(ic)+dhrwp(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp))) - txk(mxp)*dtp(ic)
           alit(2*ml,mp) = -tx(mxp)*(dxvdsh(mxp)*hrs(icxp)*usx(mxp)+  &
                (dxvdsh(mxp)*hrs(ic)+dhrsh(ic)*xvds(mxp))*(1._kdp-usx(mxp))+  &
                dxvdwh(mxp)*hrw(icxp)*uwx(mxp)+  &
                (dxvdwh(mxp)*hrw(ic)+dhrwh(ic)*xvdw(mxp))*(1._kdp-uwx(mxp)))*  &
                (p(ic)-p(icxp)) - txk(mxp)*dth(ic)
           ! ... Flow
           rlit(2*ml-1) = -tx(mxp)*  &
                ((rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*(p(icxp)-p(ic))
           ! ... Heat
           rlit(2*ml) = -tx(mxp)*  &
                ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp))*  &
                (p(icxp)-p(ic)) - txk(mxp)*(tc(icxp)-tc(ic))
        END IF
        IF (i > 1) THEN
           a1 = -tx(mxm)*  &
                ((rs(ic)*usx(mxm)+rs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (rw(ic)*uwx(mxm)+rw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                ((dxvdsp(mxm)*rs(ic)+drsp(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsp(mxm)*rs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwp(mxm)*rw(ic)+drwp(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwp(mxm)*rw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm)))
           alit(2*ml-1,mp) = alit(2*ml-1,mp) + a1
           a2 = -tx(mxm)*((dxvdsh(mxm)*rs(ic)+drsh(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsh(mxm)*rs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwh(mxm)*rw(ic)+drwh(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwh(mxm)*rw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm))
           alit(2*ml-1,mp+1) = alit(2*ml-1,mp+1) + a2
           a3 = -tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                ((dxvdsp(mxm)*hrs(ic)+dhrsp(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsp(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwp(mxm)*hrw(ic)+dhrwp(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwp(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm)))  &
                - txk(mxm)*dtp(ic)
           alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
           a4 = -tx(mxm)*((dxvdsh(mxm)*hrs(ic)+dhrsh(ic)*xvds(mxm))*usx(mxm)+  &
                dxvdsh(mxm)*hrs(icxm)*(1._kdp-usx(mxm))+  &
                (dxvdwh(mxm)*hrw(ic)+dhrwh(ic)*xvdw(mxm))*uwx(mxm)+  &
                dxvdwh(mxm)*hrw(icxm)*(1._kdp-uwx(mxm)))*(p(ic)-p(icxm))  &
                - txk(mxm)*dth(ic)
           alit(2*ml,mp) = alit(2*ml,mp) + a4
           r1 = -tx(mxm)*  &
                ((rs(ic)*usx(mxm)+rs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (rw(ic)*uwx(mxm)+rw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*(p(icxm)-p(ic))
           rlit(2*ml-1) = rlit(2*ml-1) + r1
           r2 = -tx(mxm)*  &
                ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm))*  &
                (p(icxm)-p(ic)) - txk(mxm)*(tc(icxm)-tc(ic))
           rlit(2*ml) = rlit(2*ml) + r2
        END IF
        IF (j < ny) THEN
           alit(2*ml-1,mp) = alit(2*ml-1,mp) - ty(myp)*  &
                ((rs(icyp)*usy(myp)+rs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (rw(icyp)*uwy(myp)+rw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                (dyvdsp(myp)*rs(icyp)*usy(myp)+  &
                (dyvdsp(myp)*rs(ic)+drsp(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwp(myp)*rw(icyp)*uwy(myp)+  &
                (dyvdwp(myp)*rw(ic)+drwp(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp)))
           alit(2*ml-1,mp+1) = alit(2*ml-1,mp+1) - ty(myp)*  &
                (dyvdsh(myp)*rs(icyp)*usy(myp)+  &
                (dyvdsh(myp)*rs(ic)+drsh(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwh(myp)*rw(icyp)*uwy(myp)+  &
                (dyvdwh(myp)*rw(ic)+drwh(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp))
           alit(2*ml,mp-1) = alit(2*ml,mp-1) - ty(myp)*  &
                ((hrs(icyp)*usy(myp)+hrs(ic)*(1._kdp-usy(myp)))*yvds(myp)+  &
                (hrw(icyp)*uwy(myp)+hrw(ic)*(1._kdp-uwy(myp)))*yvdw(myp)+  &
                (dyvdsp(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsp(myp)*hrs(ic)+dhrsp(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwp(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwp(myp)*hrw(ic)+dhrwp(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp))) - tyk(myp)*dtp(ic)
           alit(2*ml,mp) = alit(2*ml,mp) - ty(myp)*  &
                (dyvdsh(myp)*hrs(icyp)*usy(myp)+  &
                (dyvdsh(myp)*hrs(ic)+dhrsh(ic)*yvds(myp))*(1._kdp-usy(myp))+  &
                dyvdwh(myp)*hrw(icyp)*uwy(myp)+  &
                (dyvdwh(myp)*hrw(ic)+dhrwh(ic)*yvdw(myp))*(1._kdp-uwy(myp)))*  &
                (p(ic)-p(icyp)) - tyk(myp)*dth(ic)
        END IF
        IF (j > 1) THEN
           alit(2*ml-1,mp) = alit(2*ml-1,mp) - ty(mym)*  &
                ((rs(ic)*usy(mym)+rs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (rw(ic)*uwy(mym)+rw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                ((dyvdsp(mym)*rs(ic)+drsp(ic)*yvds(mym))*usy(mym)+  &
                dyvdsp(mym)*rs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwp(mym)*rw(ic)+drwp(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwp(mym)*rw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym)))
           alit(2*ml-1,mp+1) = alit(2*ml-1,mp+1) - ty(mym)*  &
                ((dyvdsh(mym)*rs(ic)+drsh(ic)*yvds(mym))*usy(mym)+  &
                dyvdsh(mym)*rs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwh(mym)*rw(ic)+drwh(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwh(mym)*rw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym))
           alit(2*ml,mp-1) = alit(2*ml,mp-1) - ty(mym)*  &
                ((hrs(ic)*usy(mym)+hrs(icym)*(1._kdp-usy(mym)))*yvds(mym)+  &
                (hrw(ic)*uwy(mym)+hrw(icym)*(1._kdp-uwy(mym)))*yvdw(mym)+  &
                ((dyvdsp(mym)*hrs(ic)+dhrsp(ic)*yvds(mym))*usy(mym)+  &
                dyvdsp(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwp(mym)*hrw(ic)+dhrwp(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwp(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym)))  &
                - tyk(mym)*dtp(ic)
           alit(2*ml,mp) = alit(2*ml,mp) - ty(mym)*  &
                ((dyvdsh(mym)*hrs(ic)+dhrsh(ic)*yvds(mym))*usy(mym)+  &
                dyvdsh(mym)*hrs(icym)*(1._kdp-usy(mym))+  &
                (dyvdwh(mym)*hrw(ic)+dhrwh(ic)*yvdw(mym))*uwy(mym)+  &
                dyvdwh(mym)*hrw(icym)*(1._kdp-uwy(mym)))*(p(ic)-p(icym))  &
                - tyk(mym)*dth(ic)
        END IF
        IF (k < nz) THEN
           a1 = -tz(mzp)*  &
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
           alit(2*ml-1,mp) = alit(2*ml-1,mp) + a1
           a2 = -tz(mzp)*((dzvdsh(mzp)*rs(iczp)*usz(mzp)+  &
                (dzvdsh(mzp)*rs(ic)+drsh(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwh(mzp)*rw(iczp)*uwz(mzp)+  &
                (dzvdwh(mzp)*rw(ic)+drwh(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsh(mzp)*rs(iczp)*usz(mzp)+  &
                (dzvdgsh(mzp)*rs(ic)+drsh(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwh(mzp)*rw(iczp)*uwz(mzp)+  &
                (dzvdgwh(mzp)*rw(ic)+drwh(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k))
           alit(2*ml-1,mp+1) = alit(2*ml-1,mp+1) + a2
           a3 = -tz(mzp)*  &
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
           alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
           a4 = -tz(mzp)*((dzvdsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdsh(mzp)*hrs(ic)+dhrsh(ic)*zvds(mzp))*(1._kdp-usz(mzp))+  &
                dzvdwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdwh(mzp)*hrw(ic)+dhrwh(ic)*zvdw(mzp))*(1._kdp-uwz(mzp)))*  &
                (p(ic)-p(iczp))- (dzvdgsh(mzp)*hrs(iczp)*usz(mzp)+  &
                (dzvdgsh(mzp)*hrs(ic)+dhrsh(ic)*zvdgs(mzp))*(1._kdp-usz(mzp))+  &
                dzvdgwh(mzp)*hrw(iczp)*uwz(mzp)+  &
                (dzvdgwh(mzp)*hrw(ic)+dhrwh(ic)*zvdgw(mzp))*(1._kdp-uwz(mzp)))*  &
                zzz(k)) - tzk(mzp)*dth(ic)
           alit(2*ml,mp) = alit(2*ml,mp) + a4
           r1 = -tz(mzp)*  &
                (((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
                (p(iczp)-p(ic))+  &
                ((rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
                (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))
           rlit(2*ml-1) = rlit(2*ml-1) + r1
           r2 = -tz(mzp)*  &
                (((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvds(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdw(mzp))*  &
                (p(iczp)-p(ic))+  &
                ((hrs(iczp)*usz(mzp)+hrs(ic)*(1._kdp-usz(mzp)))*zvdgs(mzp)+  &
                (hrw(iczp)*uwz(mzp)+hrw(ic)*(1._kdp-uwz(mzp)))*zvdgw(mzp))*zzz(k))  &
                - tzk(mzp)*(tc(iczp)-tc(ic))
           rlit(2*ml) = rlit(2*ml) + r2
        END IF
        IF (k > 1) THEN
           a1 = -tz(mzm)*  &
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
           alit(2*ml-1,mp) = alit(2*ml-1,mp) + a1
           a2 = -tz(mzm)*(((dzvdsh(mzm)*rs(ic)+drsh(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsh(mzm)*rs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwh(mzm)*rw(ic)+drwh(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwh(mzm)*rw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsh(mzm)*rs(ic)+drsh(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsh(mzm)*rs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwh(mzm)*rw(ic)+drwh(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwh(mzm)*rw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))
           alit(2*ml-1,mp+1) = alit(2*ml-1,mp+1) + a2
           a3 = -tz(mzm)*  &
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
           alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
           a4 = -tz(mzm)*(((dzvdsh(mzm)*hrs(ic)+dhrsh(ic)*zvds(mzm))*usz(mzm)+  &
                dzvdsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdwh(mzm)*hrw(ic)+dhrwh(ic)*zvdw(mzm))*uwz(mzm)+  &
                dzvdwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*(p(ic)-p(iczm))+  &
                ((dzvdgsh(mzm)*hrs(ic)+dhrsh(ic)*zvdgs(mzm))*usz(mzm)+  &
                dzvdgsh(mzm)*hrs(iczm)*(1._kdp-usz(mzm))+  &
                (dzvdgwh(mzm)*hrw(ic)+dhrwh(ic)*zvdgw(mzm))*uwz(mzm)+  &
                dzvdgwh(mzm)*hrw(iczm)*(1._kdp-uwz(mzm)))*zzz(k-1))  &
                - tzk(mzm)*dth(ic)
           alit(2*ml,mp) = alit(2*ml,mp) + a4
           r1 = -tz(mzm)*  &
                (((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
                (p(iczm)-p(ic))-  &
                ((rs(ic)*usz(mzm)+rs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
                (rw(ic)*uwz(mzm)+rw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))
           rlit(2*ml-1) = rlit(2*ml-1) + r1
           r2 = -tz(mzm)*  &
                (((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvds(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdw(mzm))*  &
                (p(iczm)-p(ic))-  &
                ((hrs(ic)*usz(mzm)+hrs(iczm)*(1._kdp-usz(mzm)))*zvdgs(mzm)+  &
                (hrw(ic)*uwz(mzm)+hrw(iczm)*(1._kdp-uwz(mzm)))*zvdgw(mzm))*zzz(k-1))  &
                - tzk(mzm)*(tc(iczm)-tc(ic))
           rlit(2*ml) = rlit(2*ml) + r2
        END IF
        vol = dx(i)*dy(j)*dz(k)
        IF (irad) vol = (rsq(i+1)-rsq(i))*pi*dz(k)
        a1 = -f(ic)*vol/delt
        a2 = -g(ic)*vol/delt
        a3 = -c(ic)*vol/delt + dqhp(ic)
        a4 = -d(ic)*vol/delt + dqhh(ic)
        alit(2*ml-1,mp) = alit(2*ml-1,mp) + a1
        alit(2*ml-1,mp+1) = alit(2*ml-1,mp+1) + a2
        alit(2*ml,mp-1) = alit(2*ml,mp-1) + a3
        alit(2*ml,mp) = alit(2*ml,mp) + a4
!!$        !... special output
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml-1, 2*ml-1+mp-(ibb+1), Alit(2*ml-1, mp)
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml-1, 2*ml-1+mp+1-(ibb+1), Alit(2*ml-1, mp+1)
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml, 2*ml+mp-1-(ibb+1), Alit(2*ml, mp-1)
!!$                  nnn=nnn+1
!!$                  write(60,*) 2*ml, 2*ml+mp-(ibb+1), Alit(2*ml, mp)
!!$        !...
        ! ... Add in the capacitance and source/sink terms
        r1 = (xm(ic)-xmoldt(ic))/delt - q(ic)
        r2 = (en(ic)-enoldt(ic))/delt - qh(ic)
        rlit(2*ml-1) = rlit(2*ml-1) + r1
        rlit(2*ml) = rlit(2*ml) + r2
        !    ******* Off Diagonal Block Terms  ***************
        ! ... derivatives of residual flow and transport equations w.r.t. P & H at l+1
        IF (i < nx) THEN
           nn = np(ikxp,j)
           IF (nn > 0) THEN
              nc = mp + (nn-ml)*2
              alit(2*ml-1,nc) = tx(mxp)*  &
                   ((rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                   (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                   ((dxvdsp(mxp)*rs(icxp)+drsp(icxp)*xvds(mxp))*usx(mxp)+  &
                   dxvdsp(mxp)*rs(ic)*(1._kdp-usx(mxp))+  &
                   (dxvdwp(mxp)*rw(icxp)+drwp(icxp)*xvdw(mxp))*uwx(mxp)+  &
                   dxvdwp(mxp)*rw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic)))
              alit(2*ml-1,nc+1) = tx(mxp)*  &
                   ((dxvdsh(mxp)*rs(icxp)+drsh(icxp)*xvds(mxp))*usx(mxp)+  &
                   dxvdsh(mxp)*rs(ic)*(1._kdp-usx(mxp))+  &
                   (dxvdwh(mxp)*rw(icxp)+drwh(icxp)*xvdw(mxp))*uwx(mxp)+  &
                   dxvdwh(mxp)*rw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic))
              alit(2*ml,nc-1) = tx(mxp)*  &
                   ((hrs(icxp)*usx(mxp)+hrs(ic)*(1._kdp-usx(mxp)))*xvds(mxp)+  &
                   (hrw(icxp)*uwx(mxp)+hrw(ic)*(1._kdp-uwx(mxp)))*xvdw(mxp)+  &
                   ((dxvdsp(mxp)*hrs(icxp)+dhrsp(icxp)*xvds(mxp))*usx(mxp)+  &
                   dxvdsp(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                   (dxvdwp(mxp)*hrw(icxp)+dhrwp(icxp)*xvdw(mxp))*uwx(mxp)+  &
                   dxvdwp(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic)))  &
                   + txk(mxp)*dtp(icxp)
              alit(2*ml,nc) = tx(mxp)*  &
                   ((dxvdsh(mxp)*hrs(icxp)+dhrsh(icxp)*xvds(mxp))*usx(mxp)+  &
                   dxvdsh(mxp)*hrs(ic)*(1._kdp-usx(mxp))+  &
                   (dxvdwh(mxp)*hrw(icxp)+dhrwh(icxp)*xvdw(mxp))*uwx(mxp)+  &
                   dxvdwh(mxp)*hrw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic))  &
                   + txk(mxp)*dth(icxp)
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1,2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
        IF (i > 1) THEN
           nn = np(ikxm,j)
           IF (nn > 0) THEN
              nc = mp + (nn-ml)*2
              alit(2*ml-1,nc) = tx(mxm)*  &
                   ((rs(ic)*usx(mxm)+rs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                   (rw(ic)*uwx(mxm)+rw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                   (dxvdsp(mxm)*rs(ic)*usx(mxm)+  &
                   (dxvdsp(mxm)*rs(icxm)+drsp(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                   dxvdwp(mxm)*rw(ic)*uwx(mxm)+  &
                   (dxvdwp(mxm)*rw(icxm)+drwp(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                   (p(icxm)-p(ic)))
              alit(2*ml-1,nc+1) = tx(mxm)*(dxvdsh(mxm)*rs(ic)*usx(mxm)+  &
                   (dxvdsh(mxm)*rs(icxm)+drsh(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                   dxvdwh(mxm)*rw(ic)*uwx(mxm)+  &
                   (dxvdwh(mxm)*rw(icxm)+drwh(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                   (p(icxm)-p(ic))
              alit(2*ml,nc-1) = tx(mxm)*  &
                   ((hrs(ic)*usx(mxm)+hrs(icxm)*(1._kdp-usx(mxm)))*xvds(mxm)+  &
                   (hrw(ic)*uwx(mxm)+hrw(icxm)*(1._kdp-uwx(mxm)))*xvdw(mxm)+  &
                   (dxvdsp(mxm)*hrs(ic)*usx(mxm)+  &
                   (dxvdsp(mxm)*hrs(icxm)+dhrsp(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                   dxvdwp(mxm)*hrw(ic)*uwx(mxm)+  &
                   (dxvdwp(mxm)*hrw(icxm)+dhrwp(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                   (p(icxm)-p(ic))) + txk(mxm)*dtp(icxm)
              alit(2*ml,nc) = tx(mxm)*(dxvdsh(mxm)*hrs(ic)*usx(mxm)+  &
                   (dxvdsh(mxm)*hrs(icxm)+dhrsh(icxm)*xvds(mxm))*(1._kdp-usx(mxm))+  &
                   dxvdwh(mxm)*hrw(ic)*uwx(mxm)+  &
                   (dxvdwh(mxm)*hrw(icxm)+dhrwh(icxm)*xvdw(mxm))*(1._kdp-uwx(mxm)))*  &
                   (p(icxm)-p(ic)) + txk(mxm)*dth(icxm)
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
        IF (k < nz) THEN
           nn = np(ikzp,j)
           nc = mp + (nn-ml)*2
           IF (nn > 0) THEN
              alit(2*ml-1,nc) = tz(mzp)*  &
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
              alit(2*ml-1,nc+1) = tz(mzp)*  &
                   (((dzvdsh(mzp)*rs(iczp)+drsh(iczp)*zvds(mzp))*usz(mzp)+  &
                   dzvdsh(mzp)*rs(ic)*(1._kdp-usz(mzp))+  &
                   (dzvdwh(mzp)*rw(iczp)+drwh(iczp)*zvdw(mzp))*uwz(mzp)+  &
                   dzvdwh(mzp)*rw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                   ((dzvdgsh(mzp)*rs(iczp)+drsh(iczp)*zvdgs(mzp))*usz(mzp)+  &
                   dzvdgsh(mzp)*rs(ic)*(1._kdp-usz(mzp))+  &
                   (dzvdgwh(mzp)*rw(iczp)+drwh(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                   dzvdgwh(mzp)*rw(ic)*(1._kdp-uwz(mzp)))*zzz(k))
              alit(2*ml,nc-1) = tz(mzp)*  &
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
              alit(2*ml,nc) = tz(mzp)*  &
                   (((dzvdsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvds(mzp))*usz(mzp)+  &
                   dzvdsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                   (dzvdwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdw(mzp))*uwz(mzp)+  &
                   dzvdwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*(p(iczp)-p(ic))+  &
                   ((dzvdgsh(mzp)*hrs(iczp)+dhrsh(iczp)*zvdgs(mzp))*usz(mzp)+  &
                   dzvdgsh(mzp)*hrs(ic)*(1._kdp-usz(mzp))+  &
                   (dzvdgwh(mzp)*hrw(iczp)+dhrwh(iczp)*zvdgw(mzp))*uwz(mzp)+  &
                   dzvdgwh(mzp)*hrw(ic)*(1._kdp-uwz(mzp)))*zzz(k))  &
                   + tzk(mzp)*dth(iczp)
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
        IF (k > 1) THEN
           nn = np(ikzm,j)
           nc = mp + (nn-ml)*2
           IF (nn > 0) THEN
              alit(2*ml-1,nc) = tz(mzm)*  &
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
              alit(2*ml-1,nc+1) = tz(mzm)*((dzvdsh(mzm)*rs(ic)*usz(mzm)+  &
                   (dzvdsh(mzm)*rs(iczm)+drsh(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                   dzvdwh(mzm)*rw(ic)*uwz(mzm)+  &
                   (dzvdwh(mzm)*rw(iczm)+drwh(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                   (p(iczm)-p(ic))- (dzvdgsh(mzm)*rs(ic)*usz(mzm)+  &
                   (dzvdgsh(mzm)*rs(iczm)+drsh(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                   dzvdgwh(mzm)*rw(ic)*uwz(mzm)+  &
                   (dzvdgwh(mzm)*rw(iczm)+drwh(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                   zzz(k-1))
              alit(2*ml,nc-1) = tz(mzm)*  &
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
              alit(2*ml,nc) = tz(mzm)*((dzvdsh(mzm)*hrs(ic)*usz(mzm)+  &
                   (dzvdsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvds(mzm))*(1._kdp-usz(mzm))+  &
                   dzvdwh(mzm)*hrw(ic)*uwz(mzm)+  &
                   (dzvdwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdw(mzm))*(1._kdp-uwz(mzm)))*  &
                   (p(iczm)-p(ic))- (dzvdgsh(mzm)*hrs(ic)*usz(mzm)+  &
                   (dzvdgsh(mzm)*hrs(iczm)+dhrsh(iczm)*zvdgs(mzm))*(1._kdp-usz(mzm))+  &
                   dzvdgwh(mzm)*hrw(ic)*uwz(mzm)+  &
                   (dzvdgwh(mzm)*hrw(iczm)+dhrwh(iczm)*zvdgw(mzm))*(1._kdp-uwz(mzm)))*  &
                   zzz(k-1)) + tzk(mzm)*dth(iczm)
!!$              !...
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc-(ibb+1), Alit(2*ml-1, nc)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml-1, 2*ml-1+nc+1-(ibb+1), Alit(2*ml-1, nc+1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-1-(ibb+1), Alit(2*ml, nc-1)
!!$                        nnn=nnn+1
!!$                        write(60,*) 2*ml, 2*ml+nc-(ibb+1), Alit(2*ml, nc)
!!$              !...
           END IF
        END IF
     END IF
  END DO
!!$  !...
!!$            write(60,'(1pg20.8)')  (Rlit(i),i=1,2*nby)
!!$            write(60,*) nnn, 2*nby, ibb
!!$        close(60)
!!$!*****end special output
!----Perform scaling only for 2D simulations. Scaling for 3D simulations not implemented
  IF (ny > 1) RETURN
  ! ... Scale the equations, both row and column
  ! ...     algorithm from Lapack
  !*     Get machine constants.
  smalnum = TINY(1._kdp)
  bignum = HUGE(1._kdp)
  info = 0
  !*     Compute row scale factors.
  ku = ibb
  kl = ibb
  nn = 2*nby
  ALLOCATE (r(nn),  &
       stat=a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: formeq"
     ierr(199) = .TRUE.
     RETURN
  ENDIF
  r = 0._kdp
  !..*     Find the maximum element in each row.
  DO ii = 1,nn
     DO jj = MAX(1,ii-kl),MIN(ii+ku,nn)
        r(ii) = MAX(r(ii),ABS(alit(ii,ku+1+jj-ii)))
     END DO
  END DO
  !*     Find the maximum and minimum scale factors.
  rcmin = bignum
  rcmax = 0._kdp
  rcmax = MAXVAL(r)
  rcmin = MINVAL(r)
  aamax = rcmax
  IF( rcmin == 0._kdp ) THEN
     !*        Find the first zero scale factor and return an error code.
     loc_min = MINLOC(r)
     info = loc_min(1)
     ierr(90) = .true.
     RETURN
  ELSE
     !*        Invert the scale factors.
     DO  ii = 1,nn
        r(ii) = 1._kdp/MIN(MAX(r(ii),smalnum),bignum)
     END DO
     !*        Compute ROWCND = min(R(I)) / max(R(I))
     rowcond = MAX( rcmin,smalnum ) / MIN( rcmax,bignum )
  END IF
  !*     Compute column scale factors
  cscal = 0._kdp
  !*     Find the maximum element in each column,
  !*     assuming the row scaling computed above.
  DO jj = 1,nn
     DO ii = MAX(1,jj-ku),MIN(jj+kl,nn)
        cscal(jj) = MAX(cscal(jj),ABS(alit(ii,ku+1+jj-ii))*r(ii))
     END DO
  END DO
  !*     Find the maximum and minimum scale factors.
  rcmin = bignum
  rcmax = 0._kdp
  rcmin = MINVAL(cscal(1:nn))
  rcmax = MAXVAL(cscal(1:nn))
  IF(rcmin == 0._kdp) THEN
     !*        Find the first zero scale factor and return an error code.
     loc_min = MINLOC(cscal(1:nn))
     info = loc_min(1) + nn
     ierr(90) = .true.
     RETURN
  ELSE
     !*        Invert the scale factors.
     DO jj = 1,nn
        cscal(jj) = 1._kdp/MIN(MAX(cscal(jj),smalnum),bignum)
     END DO
     !*        Compute COLCND = min(CSCAL(J)) / max(CSCAL(J))
     colcond = MAX( rcmin,smalnum ) / MIN( rcmax,bignum )
  END IF
  IF( info == 0 ) THEN
     !*           Equilibrate the matrix.
     !*        Row and column scaling
     DO  jj = 1,nn
        cj = cscal(jj)
        DO ii = MAX(1,jj-ku),MIN(jj+kl,nn)
           alit(ii,ku+1+jj-ii) = cj*r(ii)*alit(ii,ku+1+jj-ii)
        END DO
     END DO
     !*     Scale the right hand side.
     DO  ii = 1,nn
        rlit(ii) = r(ii)*rlit(ii)
     END DO
  END IF
  !*****special output for solver testing
!!$  OPEN(60,file='vaht.dat')
!!$  !...
!!$  nnn = 0
!!$  DO  ml = 1,2*nby
!!$     DO  mmr = 1,mbw
!!$        IF(ABS(alit(ml,mmr)) > 0._kdp) THEN
!!$           nnn = nnn+1
!!$           ii = ml
!!$           jj = mmr+ii-(ibb+1)
!!$           WRITE(60,*) ii, jj, alit(ml,mmr)
!!$        END IF
!!$     END DO
!!$  END DO
!!$  WRITE(60,'(1pg25.15)')  (rlit(ii),ii=1,2*nby)
!!$  WRITE(60,*) nnn, 2*nby, ibb
!!$  CLOSE(60)
!!$  !*****end special output
  DEALLOCATE(r)
END SUBROUTINE formeq
