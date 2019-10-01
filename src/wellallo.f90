SUBROUTINE wellallo
  !*****Poor subroutine name for the multiple functions
  ! ... Purpose:  To compute the distribution of mass and energy gained or lost
  ! ...        by source/sink cells and the derivatives of
  ! ...        those terms w.r.t. pressure and enthalpy.  Source/sink cells
  ! ...        include injection and production wells and 
  ! ...        the conductively heated cells along the bottom of the
  ! ...        region and the precipitation flux cells at the top
  ! ...        of the region.
  ! ...        Each column of cells is examined,
  ! ...        and if there is an active well, the mass loss rate or gain rate is
  ! ...        distributed to the cells with open intervals.  
  ! ...        Source/sink, and heat terms
  ! ...         NOTE: for flow rate terms: + =injection, - =discharge 
  ! ...         Q(IC) = mass flow rate of source/sink fluid per cell (g/s)
  ! ...         QV(IC) = volumetric flow rate of source/sink fluid per cell
  ! ...         QT(IJ) = mass flow rate of source/sink fluid per column of cells (g/s)
  ! ...         QHI(IJ) = enthalpy of injection fluid
  ! ...         CDTN(IJ) - heat flux along bottom layer of cells;
  ! ...             this flux (mW/m^2) or (erg/s-cm^2) is multiplied by the
  ! ...             area (cm^2) of each cell to give a heat flow rate (erg/s)
  ! ...         Qprecip(IJ) - precipitation flux along top layer of cells;
  ! ...             this flux (cm^3/s-cm^2) is multiplied by the
  ! ...             area (cm^2) of each cell to give a volumetric flow rate QV (cm^3/s)
  ! ... This routine works properly when there is only one point source or well
  ! ...      source term per cell.
  USE machine_constants, ONLY: kdp
  USE math_constants
  USE bc
  USE parameters          !!$ *** for flint
  USE control
  USE fdeq
  USE mesh
!!$  USE parameters
  USE source
  USE variables
  IMPLICIT NONE
  REAL(KIND=kdp) :: dsppp, dum1, dum2, dum3, dum4, dum5, dum6, dum7,  &
       dum8, dum9, dwppp, hhh, hs, hsppp, hw, hwppp, ppp, qfsum, qhhhh, qhppp
  REAL(KIND=kdp) :: rshhh, rsppp, rwhhh, rwppp, swhhh, swppp, vsppp,  &
       vwppp, xmsteam, xmwater
  INTEGER :: i, ic, icb, icpt, ij, ik, ikpt, j, k, kb
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: qf
  INTEGER :: alloc_stat
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.11 $//$Date: 2006/10/03 23:43:06 $'
  ALLOCATE(qf(nz), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: wellallo'
    ierr(177) = .TRUE.
    RETURN
  END IF
  !     ------------------------------------------------------------------
  ! ... Initialize the mass and heat source arrays to zero
  q = 0._kdp
  qh = 0._kdp
  ! ...      Set QH for conductively heated cells along base of region
  ! ...      Set QV and QH for precipitation at land surface
  ! ...      Set Q and QH for well sources
  kb = 1           ! ... Only bottom layer of cells for conduction flux
  DO  j = 1, ny
     DO  i = 1, nx
        ij = (j-1)*nx + i
        ! ... Set indices to bottom plane of cells and to the set of uppermost
        ! ...     cells within the active region
        icb = (kb-1)*nx + i + (j-1)*nx*nz
        icpt = (ktop(ij)-1)*nx + i + (j-1)*nx*nz
        ikpt = (ktop(ij)-1)*nx + i
        ! ... Apply the heat conduction flux and the preciptation flux
        ! ...      after converting to flow rates
        IF (irad) THEN        ! ... cylindrical coord
           qhcondflx(icb) = cdtn(ij)*pi*(rsq(i+1)-rsq(i))
           qh(icb) = qh(icb) + qhcondflx(icb) 
           qv(ij) = qprecip(ij)*pi*(rsq(i+1)-rsq(i))
        ELSE                  ! ... Cartesian coord
           qhcondflx(icb) = cdtn(ij)*dx(i)*dy(j)
           qh(icb) = qh(icb) + qhcondflx(icb)
           qv(ij) = qprecip(ij)*dx(i)*dy(j)
        END IF
        ! ... Do not apply precipitation flux to specified pressure cells or to
        ! ...      seepage face cells that are seeping*****
!!$        IF(ibc(ikpt,j) /= 11 .AND. ibc(ikpt,j) /= 30) THEN
!!$             q(icpt) = q(icpt) + denflux(ij)*qv(ij)
!!$             qhflux(ij) = ehflux(ij)*denflux(ij)*qv(ij)
!!$             qh(icpt) = qh(icpt) + qhflux(ij)
!!$        END IF
        ! ... Apply precipitation flux to top cells in active region
        q(icpt) = q(icpt) + denflux(ij)*qv(ij)
        qhflux(ij) = ehflux(ij)*denflux(ij)*qv(ij)
        qh(icpt) = qh(icpt) + qhflux(ij)
        ! ... No variation of source flow or heat with pressure or enthalpy for
        ! ...     heat conduction or precipitation b.c.
        dqhh(icb) = 0._kdp
        dqhp(icb) = 0._kdp
        dqhh(icpt) = 0._kdp
        dqhp(icpt) = 0._kdp
        ! ... Apply wellbore source terms
        IF(ABS(qt(ij)) > 0._kdp) THEN
           ! ... For a column of cells containing a well, calculate the mobility and
           ! ...     cell thicknesses of all open intervals and then sum in order
           ! ...     to distribute mass source or sink rates. Allocation by
           ! ...     mobility only. 
           qf = 0._kdp          ! ... Zero the layer allocation factors
           ! ... Calculate allocation factors
           DO  k = 1,nz
              ic = (k-1)*nx + i + (j-1)*nx*nz
              IF(iq(ic) == 1) THEN     ! ... A communicating cell layer
                 xmwater = rw(ic)*denw(ic)/visw(ic)
                 xmsteam = rs(ic)*dens(ic)/viss(ic)
                 ! ... Calculate allocation factor from the cell thickness and
                 ! ...     mobility. For Cartesian coordinates use a weighted average
                 ! ...     of the X and Y permeabilities
                 IF (irad) THEN
                    qf(k) = xk(ic)*dz(k)*(xmsteam+xmwater)
                 ELSE
                    qf(k) = (dx(i)*xk(ic)+dy(j)*yk(ic))/(dx(i)+dy(j))  &
                         *dz(k)*(xmsteam+xmwater)
                 END IF
              END IF
           END DO
           ! ... Sum the allocation factors for the current column of cells
           qfsum = 0._kdp
           DO  k = 1,nz
              qfsum = qfsum + qf(k)
           END DO
           ! ... Allocate mass flow rate (Q) and heat source (QH) to appropriate cells
           IF (qt(ij) > 0._kdp) THEN        ! ... Injection Well
              DO k = 1,nz
                 ic = (k-1)*nx + i + (j-1)*nx*nz
                 IF (iq(ic) == 1) THEN      ! ... A communicating cell layer
                    qwelsrc(ic) = qf(k)*qt(ij)/qfsum
                    q(ic) = q(ic) + qwelsrc(ic)
                    qhwelsrc(ic) = q(ic)*qhi(ij)
                    qh(ic) = qh(ic) + qhwelsrc(ic)
                    ! ...  Source terms are explicit
                    dqhh(ic) = 0._kdp
                    dqhp(ic) = 0._kdp
                 END IF
              END DO
           ELSE IF (qt(ij) < 0._kdp) THEN     ! ... Production Well
              DO k = 1,nz
                 ic = (k-1)*nx + i + (j-1)*nx*nz
                 IF (iq(ic) == 1) THEN     ! ... A communicating cell layer
                    ! ... Cases of single or multi-phase production
                    IF (ind(ic) /= 2) THEN
                       ! ... Compressed water or saturated steam or air-water
                       qwelsrc(ic) = qf(k)*qt(ij)/qfsum
                       q(ic) = q(ic) + qwelsrc(ic)
                       qhwelsrc(ic) = q(ic)*h(ic)
                       qh(ic) = qh(ic) + qhwelsrc(ic)
                       dqhh(ic) = q(ic)     ! ... Heat sink depends on enthalpy
                       dqhp(ic) = 0._kdp
                    ELSE
                       ! ... Two-Phase Region
                       qwelsrc(ic) = qf(k)*qt(ij)/qfsum
                       q(ic) = q(ic) + qwelsrc(ic)
                       ! ... Calculate enthalpies of saturated water (HW) and steam (HS) at 
                       ! ...     the cell pressure
                       CALL phsbndy1(p(ic),hw,hs)
                       ! ... Apply heat sink in proportion to mobility of each phase
                       qhwelsrc(ic) = q(ic)*  &
                            (hw*rw(ic)*denw(ic)*viss(ic) + hs*rs(ic)*dens(ic)*visw(ic))  &
                            /(rw(ic)*denw(ic)*viss(ic) + rs(ic)*dens(ic)*visw(ic))
                       qh(ic) = qh(ic) + qhwelsrc(ic)
                       ! ... Calculate the derivative of QH with respect to enthalpy
                       ! ...      Use a finite perturbation approximation. Change enthalpy by 
                       ! ...      .01% and use new heat source value to determine derivatives;
                       ! ...      The enthalpy perturbation is made so as to avoid crossing the 
                       ! ...      two-phase curve,
                       ! ...    i.e. increase enthalpy by 1.0001X unless water sat'n is < 0.1,
                       ! ...    then decrease enthalpy by 0.9999X
                       hhh = h(ic)*1.0000001_kdp
                       IF (satnw(ic) < 0.1_kdp) hhh = h(ic)*0.999999_kdp
                       ! ... It is not necessary to recalculate enthalpies, densities or viscosities
                       ! ... of saturated water and steam, or temperature since they do not vary
                       ! ... with enthalpy in the 2 phase zone (i.e. derivatives of density,
                       ! ... viscosity, and enthalpy of steam and water, & temp w.r.t. enthalpy= 0)
                       ! ... Recalculate water saturation (SWHHH) using perturbed enthalpy (HHH)
                       swhhh = dens(ic)*(hs-hhh)/(dens(ic)*(hs-hhh)+denw(ic)*(hhh-hw))
                       ! ... Recalculate relative permeabilities for water (RWHHH) and steam (RSHHH)
                       ! ...    using new water saturation (SWN) based on perturbed enthalpy (HHH)
                       CALL relperm(kodrp, p(ic), swhhh, 0._kdp, 0._kdp,  &
                            rwhhh, dum1, dum2, rshhh, dum3, dum4, ind(ic))
                       qhhhh = q(ic)*  &
                            (hw*rwhhh*denw(ic)*viss(ic) + hs*rshhh*dens(ic)*visw(ic))  &
                            /(rwhhh*denw(ic)*viss(ic) + rshhh*dens(ic)*visw(ic))
                       ! ... If node is on bottom row, include the conductive heat b.c.
                       IF (k == 1) THEN
                          IF (irad) THEN
                             qhhhh = qhhhh + cdtn(ij)*pi*(rsq(i+1)-rsq(i))
                          ELSE
                             qhhhh = qhhhh + cdtn(ij)*dx(i)*dy(j)
                          END IF
                       END IF
                       ! ... If node is at uppermost active cell, include the advective heat  
                       ! ...      flow from precipitation
                       IF (k == ikpt) qhhhh = qhhhh + qhflux(ij)
                       ! ... Finite perturbation approximation for derivative with respect to
                       ! ...      enthalpy
                       dqhh(ic) = (qhhhh-qh(ic))/(hhh-h(ic))
                       ! ... Calculate the derivative of QH with respect to pressure
                       ! ...      Use a finite perturbation approximation. Change pressure by 
                       ! ...      .01% and use new heat source value to determine derivatives;
                       ! ...      The pressure perturbation is made so as to avoid crossing the 
                       ! ...      two-phase curve,
                       ! ...      i.e. if water sat'n is < 0.1 and pressure is < 70 bars, 
                       ! ...      then increase the pressure
                       ppp = p(ic)*0.9999999_kdp
                       IF (satnw(ic) < 0.1_kdp .AND. p(ic) < 7.0e7_kdp) ppp = p(ic)*1.0000001_kdp
                       IF (p(ic) < 0.6e6_kdp) ppp = p(ic)*1.0000001_kdp
                       ! ... Calculate the enthalpy of sat'd water and steam and
                       ! ...      the density, temperature, and viscosity of
                       ! ...      sat'd water and sat'd steam for the two-phase region
                       CALL phsbndy3(ppp, hwppp, hsppp, dum1, dum2,  &
                            dwppp, dum3, vwppp, dum4, dsppp, dum5, vsppp, dum6,  &
                            dum7, dum8, dum9)
                       ! ... Recalculate water saturation (SWPPP) using modified pressure (PPP)
                       swppp = dsppp*(hsppp-h(ic))/(dsppp*(hsppp-h(ic))+dwppp*(h(ic)-hwppp))
                       ! ... Recalculate relative permeabilities for water (RWPPP) and steam (RSPPP)
                       ! ...    using new water satn (SWPPP) based on perturbed pressure (PPP)
                       CALL relperm(kodrp, p(ic), swppp, 0._kdp, 0._kdp,  &
                            rwppp, dum1, dum2, rsppp, dum3, dum4, ind(ic))
                       qhppp = q(ic)*  &
                            (hwppp*rwppp*dwppp*vsppp+hsppp*rsppp*dsppp*vwppp)  &
                            /(rwppp*dwppp*vsppp+rsppp*dsppp*vwppp)
                       ! ... If node is on bottom row, include the conductive heat b.c.
                       IF (k == 1) THEN
                          IF (irad) THEN
                             qhppp = qhppp + cdtn(ij)*pi*(rsq(i+1)-rsq(i))
                          ELSE
                             qhppp = qhppp + cdtn(ij)*dx(i)*dy(j)
                          END IF
                       END IF
                       ! ... If node is at uppermost active cell, include the advective heat  
                       ! ...      flow from precipitation
                       IF (k == ikpt) qhppp = qhppp + qhflux(ij)
                       ! ... Finite perturbation approximation for derivative with respect to
                       ! ...      pressure
                       dqhp(ic) = (qhppp-qh(ic))/(ppp-p(ic))
                    END IF
                 END IF
              END DO
           END IF
        END IF
     END DO
  END DO
  ! ... Now apply the nodal point source mass and heat flow rates
  DO  ic=1,nxyz
     call ictoijk(ic,i,j,k,nx,nz)
     ij = (j-1)*nx + i
     IF(iq(ic) == 1) THEN        ! ... For a node that has a source, allocate it
        IF (qnodsrc(ic) > 0._kdp) THEN          ! ... Source Node
           q(ic) = q(ic) + qnodsrc(ic)
           qhnodsrc(ic) = qnodsrc(ic)*qhinodsrc(ic)
           qh(ic) = qh(ic) + qhnodsrc(ic)
           ! ... Source terms are explicit
           dqhh(ic) = 0._kdp
           dqhp(ic) = 0._kdp
        ELSE IF (qnodsrc(ic) < 0._kdp) THEN     ! ... Sink Node
           ! ... Cases of single or multi-phase production
           IF (ind(ic) /= 2) THEN
              ! ... Compressed water or saturated steam or air-water
              q(ic) = q(ic) + qnodsrc(ic)
              qhnodsrc(ic) = qnodsrc(ic)*h(ic)
              qh(ic) = qh(ic) + qhnodsrc(ic)
              dqhh(ic) = q(ic)
              dqhp(ic) = 0._kdp
           ELSE
              ! ... Two Phase Region 
              q(ic) = q(ic) + qnodsrc(ic)
              ! ... Calculate enthalpies of saturated water (HW) and steam (HS) at 
              ! ...       the cell pressure
              CALL phsbndy1(p(ic),hw,hs)
              ! ... Apply heat sink in proportion to mobility of each phase
              qhnodsrc(ic) = qnodsrc(ic)*  &
                   (hw*rw(ic)*denw(ic)*viss(ic) + hs*rs(ic)*dens(ic)*visw(ic))  &
                   /(rw(ic)*denw(ic)*viss(ic) + rs(ic)*dens(ic)*visw(ic))
              qh(ic) = qh(ic) + qhnodsrc(ic)
              ! ... Calculate the derivative of QH with respect to enthalpy
              ! ...      Use a finite perturbation approximation. Change enthalpy by 
              ! ...      .01% and use new heat source value to determine derivatives;
              ! ...      The enthalpy perturbation is made so as to avoid crossing the 
              ! ...      two-phase curve,
              ! ...    i.e. increase enthalpy by 1.0001X unless water sat'n is < 0.1,
              ! ...    then decrease enthalpy by 0.9999X
              hhh = h(ic)*1.0000001_kdp
              IF (satnw(ic) < 0.1_kdp) hhh = h(ic)*0.999999_kdp
              ! ... It is not necessary to recalculate enthalpies, densities or viscosities
              ! ... of saturated water and steam, or temperature since they do not vary
              ! ... with enthalpy in the 2 phase field (i.e. derivatives of density,
              ! ... viscosity, and enthalpy of steam and water, & temp w.r.t. enthalpy= 0)
              ! ... recalculate water saturation (SWHHH) using modified enthalpy (HHH)
              swhhh = dens(ic)*(hs-hhh)/(dens(ic)*(hs-hhh)+denw(ic)*(hhh-hw))
              ! ... Recalculate relative permeabilities for water (RWHHH) and steam (RSHHH)
              ! ...    using new water saturation (SWN) based on modified enthalpy (HHH)
              CALL relperm(kodrp, p(ic), swhhh, 0._kdp, 0._kdp,  &
                   rwhhh, dum1, dum2, rshhh, dum3, dum4, ind(ic))
              qhhhh = qnodsrc(ic)*  &
                   (hw*rwhhh*denw(ic)*viss(ic)+hs*rshhh*dens(ic)*visw(ic))  &
                   /(rwhhh*denw(ic)*viss(ic)+rshhh*dens(ic)*visw(ic))
              ! ... If node is on bottom row, include the conductive heat b.c.
              IF (k == 1) THEN
                 IF (irad) THEN
                    qhhhh = qhhhh + cdtn(ij)*pi*(rsq(i+1)-rsq(i))
                 ELSE
                    qhhhh = qhhhh + cdtn(ij)*dx(i)*dy(j)
                 END IF
              END IF
              ! ... If node is at uppermost active cell, include the advective heat  
              ! ...      flow from precipitation
              IF (k == ikpt) qhhhh = qhhhh + qhflux(ij)
              ! ... Finite perturbation approximation for derivative with respect to
              ! ...      enthalpy
              dqhh(ic) = (qhhhh-qh(ic))/(hhh-h(ic))
              ! ... Calculate the derivative of QH with respect to pressure
              ! ...      Use a finite perturbation approximation. Change pressure by 
              ! ...      .01% and use new heat source value to determine derivatives;
              ! ...      The pressure perturbation is made so as to avoid crossing the 
              ! ...      two-phase curve,
              ! ...      i.e. if water sat'n is < 0.1 and pressure is < 70 bars, 
              ! ...      then increase the pressure
              ppp = p(ic)*0.9999999_kdp
              IF (satnw(ic) < 0.1_kdp .AND. p(ic) < 7.0e7_kdp) ppp = p(ic)*1.0000001_kdp
              IF (p(ic) < 0.6e6_kdp) ppp = p(ic)*1.0000001_kdp
              ! ... Calculate the enthalpy of sat'd water and steam and
              ! ...      the density, temperature, and viscosity of
              ! ...      sat'd water and sat'd steam for the two-phase region
              CALL phsbndy3(ppp, hwppp, hsppp, dum1, dum2,  &
                   dwppp, dum3, vwppp, dum4, dsppp, dum5, vsppp, dum6,  &
                   dum7, dum8, dum9)
              ! ... Recalculate water saturation (SWPPP) using modified pressure (PPP)
              swppp = dsppp*(hsppp-h(ic))/  &
                   (dsppp*(hsppp-h(ic))+dwppp*(h(ic)-hwppp))
              ! ... Recalculate relative permeabilities for water (RWPPP) and steam (RSPPP)
              ! ...      using new water satn (SWPPP) based on modified pressure (PPP)
              CALL relperm(kodrp, p(ic), swppp, 0._kdp, 0._kdp,  &
                   rwppp, dum1, dum2, rsppp, dum3, dum4, ind(ic))
              qhppp = q(ic)*  &
                   (hwppp*rwppp*dwppp*vsppp+hsppp*rsppp*dsppp*vwppp)  &
                   /(rwppp*dwppp*vsppp+rsppp*dsppp*vwppp)
              ! ... If node is on bottom row, include the conductive heat b.c.
              IF (k == 1) THEN
                 IF (irad) THEN
                    qhppp = qhppp + cdtn(ij)*pi*(rsq(i+1)-rsq(i))
                 ELSE
                    qhppp = qhppp + cdtn(ij)*dx(i)*dy(j)
                 END IF
              END IF
              ! ... If node is at uppermost active cell, include the advective heat  
              ! ...      flow from precipitation
              IF (k == ikpt) qhppp = qhppp + qhflux(ij)
              ! ... Finite perturbation approximation for derivative with respect to
              ! ...      pressure
              dqhp(ic) = (qhppp-qh(ic))/(ppp-p(ic))
           END IF
        END IF
     END IF
  END DO
  DEALLOCATE(qf)
  ! ... *** The following is now done in subroutine formeq or asem_nr_csr
!!$  ! ... Now apply the seepage surface sink mass and heat flow rates
!!$  IF(nseep > 0) THEN
!!$     DO ic=1,nxyz
!!$        call ictoijk(ic,i,j,k,nx,nz)
!!$        ij = (j-1)*nx + i
!!$        ik = (k-1)*nx + i
!!$        IF(ibc(ik,j) == 30 .AND. seeping(ik,j)) THEN
!!$           ! ... For a seeping node, calculate the heat sink term.
!!$           ! ... Only the case of single-phase seepage is allowed.
!!$           ! ...      Compressed water or air-water
!!$!...              q(ic) = q(ic) + qnodsrc(ic)  ! ... No source for sp. pressure seeping cell
!!$              qh(ic) = qh(ic) + qseep(ic)*h(ic)
!!$              dqhh(ic) = qseep(ic)
!!$              dqhp(ic) = 0._kdp
!!$        END IF
!!$     END DO
!!$  END IF
END SUBROUTINE wellallo
