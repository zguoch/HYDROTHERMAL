SUBROUTINE prpty(ilfast)
  ! ... Purpose:  To compute the physical properties of the fluid:
  ! ...           density, viscosity, relative permeability, saturation,
  ! ...           and phase region for both the water and steam phases.  In addition,
  ! ...           the derivatives of these variables are computed with respect to
  ! ...           enthalpy and pressure.  The variables and derivatives are calculated at
  ! ...           the nodes.  The arithmetic average of the
  ! ...           density and viscosity are determined for the cell faces.
  ! ...           The derivatives of these cell face variables are
  ! ...           incomplete because they lack the weighting
  ! ...           of the averaged variable; this weighting is inserted
  ! ...           into the Newton equation during assembly.
  USE machine_constants, ONLY: kdp
  USE f_units
  USE parameters          !!$*** for flint
  USE control
  USE i_c
  USE fdeq
  USE mesh
!!$  USE parameters
  USE variables
  USE units
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ilfast
  ! ...   ILFAST (input) = true(1)  - set properties for only the active nodes
  ! ...                  = false(0) - set properties for all nodes
  !
  INTEGER :: i, ic, ik, j, k
  LOGICAL :: errflag, lfast
  REAL(kind=kdp) :: dum1, hs, hw, pdsh, pdsp, pdwh, pdwp, phsp, phwp,  &
       pswh, pswp, pvsh, pvsp, pvwh, pvwp, usw, vsw
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4.1.3 $//$Date: 2006/10/03 23:24:37 $'
  !     ------------------------------------------------------------------
  lfast = ilfast == 1
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF (lfast .AND. np(ik,j) <= 0) CYCLE       ! ... skip non-active nodes when appropriate
     ! ...         if node not in domain set density, viscosity, and
     ! ...              saturation terms for later arithmetic averaging
     ! ...              then skip to end of loop
     IF (np(ik,j) == 0) THEN
        denw(ic) = 0._kdp
        dens(ic) = 0._kdp
        visw(ic) = 1.0e20_kdp
        viss(ic) = 1.0e20_kdp
        satnw(ic) = 1._kdp
        rw(ic) = 0._kdp
        rs(ic) = 0._kdp
        CYCLE
     END IF
     IF(ind(ic) /= 5) THEN
        ! ... Calculate new porosity; but not for air-water cells
        phi(ic) = phiinit(ic)*(1._kdp + beta(ic)*(p(ic) - pinit(ic)))
        IF (phi(ic) > 1._kdp .OR. phi(ic) < 0._kdp) THEN
           ! ... Error stop if porosity is out of range
           WRITE (fustdout, 9005) i, j, k, ic, phi(ic)
           WRITE (fuclog, 9005) i, j, k, ic, phi(ic)
9005       FORMAT (/  &
                ' ***  Porosity, PHI, out of range at I, J, K  ',  &
                3I3 / 10X, '  IC=', i4, '   PHI(IC)=', g13.6)
           ierr(101) = .TRUE.
           !..              return
        END IF
     END IF
     ! ...     go to the appropriate thermodynamic phase region for the cell as
     ! ...          determined in PHREG and flagged in IND(IC)
     IF (ind(ic) == 1 .AND. p(ic) < 220.55e6_kdp) THEN
        ! ... ****** ****** ****** COMPRESSED WATER REGION ******* ****** ******
        ! ...       determine density, temperature, and viscosity of water
        ! ...         and the derivatives of these functions w.r.t. H and P
        hw = 0._kdp
        hs = 0._kdp
        phwp = 0._kdp
        phsp = 0._kdp
        pswp = 0._kdp
        pswh = 0._kdp
        CALL tblookup(p(ic), h(ic),denw(ic), pdwh, pdwp,  &
             tc(ic), dth(ic), dtp(ic), visw(ic), pvwh, pvwp, errflag)
        IF (p(ic) > 220.0e6_kdp) THEN
           ! ...   NOTE:  water and steam properties are made identical for pressures
           ! ...          between 220.0 and 220.55 (Pcrit) (dyne/cm^2) to smooth transition into
           ! ...          supercritical region
           dens(ic) = denw(ic)
           viss(ic) = visw(ic)
           pdsp = pdwp
           pdsh = pdwh
           pvsp = pvwp
           pvsh = pvwh
           ! ... set water & steam relative permeabilities and derivatives
           CALL relperm2(p(ic), h(ic), rw(ic), drwp(ic),drwh(ic),  &
                rs(ic), drsp(ic), drsh(ic))
        ELSE
           dens(ic) = 0._kdp
           viss(ic) = 1.0e20_kdp
           pdsp = 0._kdp
           pdsh = 0._kdp
           pvsp = 0._kdp
           pvsh = 0._kdp
           ! ... set water & steam relative permeabilities and derivatives
           rw(ic) = 1._kdp
           rs(ic) = 0._kdp
           drwp(ic) = 0._kdp
           drwh(ic) = 0._kdp
           drsp(ic) = 0._kdp
           drsh(ic) = 0._kdp
        END IF
        ! ...      set water and steam saturation
        satnw(ic) = 1._kdp
        ! ... calculate the heat times relative permeability terms
        hrw(ic) = h(ic)*rw(ic)
        hrs(ic) = h(ic)*rs(ic)
        ! ...     derivatives of heat times relative permeability terms
        ! ...          water
        dhrwh(ic) = h(ic)*drwh(ic) + rw(ic)
        dhrwp(ic) = h(ic)*drwp(ic)
        ! ...          steam
        dhrsh(ic) = h(ic)*drsh(ic) + rs(ic)
        dhrsp(ic) = h(ic)*drsp(ic)
        ! ...     source/sink heat terms
        ! ...         The variables for source and sink terms Q(IC) and QH(IC)
        ! ...         and their derivatives DQHH & DQHP are set in WELLALLO.
        ! ...     derivatives of viscosity and density terms are calculated
        ! ...       for the interblock region in subroutines avgvalue and avgprpty
        ! ... compute derivatives of capacitance in governing equations
        ! ...      Energy equation
        c(ic) = pdwp*phi(ic)*h(ic) +  &
             (1._kdp-phi(ic))*df(ic)*phfwt(ic)*dtp(ic) +  &
             phiinit(ic)*beta(ic)*denw(ic)*h(ic) -  &
             phiinit(ic)*beta(ic)*df(ic)*hrock(ic)
        d(ic) = phi(ic)*denw(ic) + pdwh*phi(ic)*h(ic) +  &
             (1._kdp-phi(ic))*df(ic)*phfwt(ic)*dth(ic) +  &
             (1._kdp-phi(ic))*df(ic)*phfwtdt(ic)*dth(ic)*(273.15_kdp+tc(ic))
        ! ...      Flow equation
        f(ic) = phi(ic)*pdwp + denw(ic)*beta(ic)*phiinit(ic)
        g(ic) = phi(ic)*pdwh
        CYCLE
     ELSE IF (ind(ic) == 2) THEN
        ! ...  ******* ****** ****** TWO-PHASE REGION ******* ****** ******
        ! ...       determine enthalpy of sat'd water and steam and
        ! ...         the density, temperature, and viscosity of
        ! ...         sat'd water and sat'd steam for the two-phase region
        ! ...         and the derivatives of these functions w.r.t. H and P
        CALL phsbndy3(p(ic), hw, hs, phwp, phsp,  &
             denw(ic), pdwp, visw(ic), pvwp, dens(ic), pdsp, viss(ic), pvsp,  &
             tc(ic), dtp(ic), dum1)
        pdwh = 0._kdp
        pdsh = 0._kdp
        pvwh = 0._kdp
        pvsh = 0._kdp
        dth(ic) = 0._kdp
        ! ...     calculate water saturation (SATNW) and its derivatives
        usw = dens(ic)*(hs-h(ic))
        vsw = dens(ic)*(hs-h(ic)) + denw(ic)*(h(ic)-hw)
        satnw(ic) = usw/vsw
        satnw(ic) = MIN(1._kdp, satnw(ic))
        satnw(ic) = MAX(0._kdp, satnw(ic))
        pswp = (pdsp*(hs-h(ic))+dens(ic)*phsp)/vsw -  &
             (pdsp*(hs-h(ic))+dens(ic)*phsp  &
             +pdwp*(h(ic)-hw)-denw(ic)*phwp)*usw/(vsw**2)
        pswh = -1._kdp*dens(ic)/vsw - (denw(ic)-dens(ic))*usw/(vsw**2)
        ! ...    calculate relative permeability for water (RW) and steam (RS)
        CALL relperm(kodrp, p(ic), satnw(ic), pswp, pswh,  &
             rw(ic), drwp(ic), drwh(ic), rs(ic), drsp(ic), drsh(ic), ind(ic))
        ! ... calculate heat times relative permeability terms
        hrw(ic) = hw*rw(ic)
        hrs(ic) = hs*rs(ic)
        ! ...     derivatives of heat * relative perm. terms
        ! ...          water
        ! ...        if the sat'n of water = 1 (this happens for two phase nodes
        ! ...          when nodes just barely cross the phase line - see phreg.f)
        ! ...          calculate DHRWH(IC) the same way its done for compressed water
        dhrwp(ic) = hw*drwp(ic) + rw(ic)*phwp
        dhrwh(ic) = hw*drwh(ic)
        IF (satnw(ic) == 1._kdp) dhrwh(ic) = 1._kdp
        ! ...          steam
        ! ...        if the sat'n of water = 0 (this happens for two phase nodes
        ! ...          when nodes just barely cross the phase line - see phreg.f)
        ! ...          calculate DHRSH(IC) the same way its done for steam
        dhrsp(ic) = hs*drsp(ic) + rs(ic)*phsp
        dhrsh(ic) = hs*drsh(ic)
        IF (satnw(ic) == 0._kdp) dhrsh(ic) = 1._kdp
        ! ...     source/sink heat terms
        ! ...         The variables for source and sink terms Q(IC) and QH(IC)
        ! ...         and their derivatives DQHH & DQHP are set in WELLALLO.
        ! ...     derivatives of viscosity and density terms are calculated
        ! ...       for the interblock region in subroutines avgvalue and avgprpty
        ! ... compute derivatives of capacitance in governing equations
        c(ic) = phiinit(ic)*beta(ic) *((1._kdp-satnw(ic))*dens(ic)*hs  &
             +satnw(ic)*denw(ic)*hw-df(ic)*phfwt(ic)*(273.15_kdp+tc(ic)))  &
             - pswp*phi(ic)*dens(ic)*hs + pdsp*phi(ic)*(1._kdp-satnw(ic))*hs  &
             + phsp*phi(ic)*(1._kdp-satnw(ic))*dens(ic)  &
             + pswp*phi(ic)*denw(ic)*hw + pdwp*phi(ic)*satnw(ic)*hw  &
             + phwp*phi(ic)*satnw(ic)*denw(ic)  &
             + (1._kdp-phi(ic))*df(ic)*phfwt(ic)*dtp(ic)
        d(ic) = pswh*phi(ic)*(denw(ic)*hw-dens(ic)*hs)
        f(ic) = phiinit(ic)*beta(ic)*(satnw(ic)*denw(ic)  &
             +(1._kdp-satnw(ic))*dens(ic)) + pswp*phi(ic)*(denw(ic)-dens(ic))  &
             + phi(ic)*(pdwp*satnw(ic)+pdsp*(1._kdp-satnw(ic)))
        g(ic) = pswh*phi(ic)*(denw(ic)-dens(ic))
        CYCLE
     ELSE IF (ind(ic) == 3) THEN
        ! ... ****** ******* ****** SUPER-HEATED STEAM REGION ****** ******* ******
        ! ...   NOTE:  water and steam properties are made identical for pressures
        ! ...          between 220.0 and 220.55 (Pcrit) to smooth changes into
        ! ...          supercritical region
        ! ...       determine density, temperature, and viscosity of steam
        ! ...         and the derivatives of these functions w.r.t. H and P
        hw = 0._kdp
        hs = 0._kdp
        phwp = 0._kdp
        phsp = 0._kdp
        pswp = 0._kdp
        pswh = 0._kdp
        CALL tblookup(p(ic), h(ic), dens(ic), pdsh, pdsp,  &
             tc(ic), dth(ic), dtp(ic), viss(ic), pvsh, pvsp, errflag)
        IF (p(ic) > 220.0e6_kdp) THEN
           denw(ic) = dens(ic)
           visw(ic) = viss(ic)
           pdwp = pdsp
           pdwh = pdsh
           pvwp = pvsp
           pvwh = pvsh
        ELSE
           denw(ic) = 0._kdp
           visw(ic) = 1.0e20_kdp
           pdwp = 0._kdp
           pdwh = 0._kdp
           pvwp = 0._kdp
           pvwh = 0._kdp
        END IF
        ! ...      set water and steam saturation
        satnw(ic) = 0._kdp
        ! ...      set water & steam relative permeabilities and derivatives
        IF (p(ic) > 220.0e6_kdp) THEN
           CALL relperm2(p(ic), h(ic), rw(ic), drwp(ic), drwh(ic),  &
                rs(ic), drsp(ic), drsh(ic))
        ELSE
           ! ...          water
           rw(ic) = 0._kdp
           drwh(ic) = 0._kdp
           drwp(ic) = 0._kdp
           ! ...          steam
           rs(ic) = 1._kdp
           drsh(ic) = 0._kdp
           drsp(ic) = 0._kdp
        END IF
        ! ...      calculated heat times relative permeability terms
        hrw(ic) = h(ic)*rw(ic)
        hrs(ic) = h(ic)*rs(ic)
        ! ...     derivatives of heat * relative perm. terms
        ! ...          water
        dhrwh(ic) = h(ic)*drwh(ic) + 1._kdp*rw(ic)
        dhrwp(ic) = h(ic)*drwp(ic)
        ! ...          steam
        dhrsh(ic) = h(ic)*drsh(ic) + 1._kdp*rs(ic)
        dhrsp(ic) = h(ic)*drsp(ic)
        ! ...     source/sink heat terms
        ! ...         The variables for source and sink terms Q(IC) and QH(IC)
        ! ...         and their derivatives DQHH & DQHP are set in WELLALLO.
        ! ...     derivatives of viscosity and density terms are calculated
        ! ...       for the interblock region in subroutines avgvalue and avgprpty
        ! ... compute derivatives of capacitance in governing equations
        c(ic) = pdsp*phi(ic)*h(ic) +  &
             (1._kdp-phi(ic))*df(ic)*phfwt(ic)*dtp(ic) +  &
             phiinit(ic)*beta(ic)*dens(ic)*h(ic) -  &
             phiinit(ic)*beta(ic)*df(ic)*hrock(ic)
        !    &   PHIINIT(IC)*BETA(IC)*DF(IC)*PHFWT(IC)*(273.15_kdp+TC(IC))
        d(ic) = phi(ic)*dens(ic) + pdsh*phi(ic)*h(ic) +  &
             (1._kdp-phi(ic))*df(ic)*phfwt(ic)*dth(ic) +  &
             (1._kdp-phi(ic))*df(ic)*phfwtdt(ic)*dth(ic)*(273.15_kdp+tc(ic))
        f(ic) = phi(ic)*pdsp + dens(ic)*beta(ic)*phiinit(ic)
        g(ic) = phi(ic)*pdsh
        CYCLE
     ELSE IF((ind(ic) == 4 .OR. ind(ic) == 1) .AND. p(ic) >= 220.55e6_kdp) THEN
        ! ... ****** ******* ****** SUPER-CRITICAL FLUID ****** ******* ******
        ! ...       determine density, temperature, and viscosity of steam
        ! ...         and the derivatives of these functions w.r.t. H and P
        hw = 0._kdp
        hs = 0._kdp
        phwp = 0._kdp
        phsp = 0._kdp
        pswp = 0._kdp
        pswh = 0._kdp
        CALL tblookup(p(ic), h(ic), dens(ic), pdsh, pdsp,  &
             tc(ic), dth(ic), dtp(ic), viss(ic), pvsh, pvsp, errflag)
        denw(ic) = dens(ic)
        visw(ic) = viss(ic)
        pdwp = pdsp
        pdwh = pdsh
        pvwp = pvsp
        pvwh = pvsh
        ! ...      set water and steam saturation
        satnw(ic) = 0.5_kdp
        ! ...      set water & steam relative permeabilities
        CALL relperm2(p(ic), h(ic), rw(ic), drwp(ic), drwh(ic),  &
             rs(ic), drsp(ic), drsh(ic))
        ! ...      calculated heat times relative permeability terms
        hrw(ic) = h(ic)*rw(ic)
        hrs(ic) = h(ic)*rs(ic)
        ! ...     derivatives of heat * relative perm. terms
        ! ...          water
        dhrwh(ic) = h(ic)*drwh(ic) + 1._kdp*rw(ic)
        dhrwp(ic) = h(ic)*drwp(ic)
        ! ...          steam
        dhrsh(ic) = h(ic)*drsh(ic) + 1._kdp*rs(ic)
        dhrsp(ic) = h(ic)*drsp(ic)
        ! ...     source/sink heat terms
        ! ...         The variables for source and sink terms Q(IC) and QH(IC)
        ! ...         and their derivatives DQHH & DQHP are set in WELLALLO.
        ! ...     derivatives of viscosity and density terms are calculated
        ! ...       for the interblock region in subroutines avgvalue and avgprpty
        ! ... compute derivatives of capacitance in governing equations
        c(ic) = pdsp*phi(ic)*h(ic) +  &
             (1._kdp-phi(ic))*df(ic)*phfwt(ic)*dtp(ic) +  &
             phiinit(ic)*beta(ic)*dens(ic)*h(ic) -  &
             phiinit(ic)*beta(ic)*df(ic)*hrock(ic)
        d(ic) = phi(ic)*dens(ic) + pdsh*phi(ic)*h(ic) +  &
             (1._kdp-phi(ic))*df(ic)*phfwt(ic)*dth(ic) +  &
             (1._kdp-phi(ic))*df(ic)*phfwtdt(ic)*dth(ic)*(273.15_kdp+tc(ic))
        f(ic) = phi(ic)*pdsp + dens(ic)*beta(ic)*phiinit(ic)
        g(ic) = phi(ic)*pdsh
     ELSE IF (ind(ic) == 5) THEN
        ! ...  ******* ****** ****** Air-Water Two-phase Region ******* ****** ******
        IF (p(ic) < 220.55e6_kdp) THEN
           ! ...       determine density, temperature, and viscosity of water
           ! ...         and the derivatives of these functions w.r.t. H and P
           hw = h(ic)
           hs = 0._kdp
           phwp = 0._kdp
           phsp = 0._kdp
           pswp = 0._kdp
           pswh = 0._kdp
           CALL tblookup(p(ic), h(ic),denw(ic), pdwh, pdwp,  &
                tc(ic), dth(ic), dtp(ic), visw(ic), pvwh, pvwp, errflag)
           ! ... Steam density and derivatives are set to zero, viscosity set to very
           ! ...      large value
           dens(ic) = 0._kdp
           viss(ic) = 1.0e20_kdp
           pdsp = 0._kdp
           pdsh = 0._kdp
           pvsp = 0._kdp
           pvsh = 0._kdp
           ! ...     calculate water saturation (SATNW) and its derivatives
           CALL satur_w(kodrp,p(ic),h(ic),satnw(ic),pswp,pswh)
           ! ...    calculate relative permeability for water (RW)
           ! ... Steam values are returned but not used, so they are reset to zero
           ssr=0._kdp
           CALL relperm(kodrp, p(ic), satnw(ic), pswp, pswh,  &
                rw(ic), drwp(ic), drwh(ic), rs(ic), drsp(ic), drsh(ic), ind(ic))
           rs(ic) = 0._kdp
           ! ... calculate the enthalpy-relative permeability product
           hrw(ic) = h(ic)*rw(ic)
           hrs(ic) = 0._kdp
           ! ...     derivatives of heat times relative permeability terms
           ! ...          water
           dhrwh(ic) = h(ic)*drwh(ic) + rw(ic)
           dhrwp(ic) = h(ic)*drwp(ic)
           ! ...          steam
           dhrsh(ic) = 0._kdp
           dhrsp(ic) = 0._kdp
           ! ...     source/sink heat terms
           ! ...         The variables for source and sink terms Q(IC) and QH(IC)
           ! ...         and their derivatives DQHH & DQHP are set in WELLALLO.
           ! ...     derivatives of viscosity and density terms are calculated
           ! ...       for the interblock region in subroutines avgvalue and avgprpty
           ! ... compute derivatives of capacitance in governing equations
           c(ic) = pswp*phi(ic)*denw(ic)*hw + pdwp*phi(ic)*satnw(ic)*hw + &
                (1._kdp-phi(ic))*df(ic)*phfwt(ic)*dtp(ic)
           d(ic) = pswh*phi(ic)*denw(ic)*hw + pdwh*phi(ic)*satnw(ic)*hw + &
                phi(ic)*satnw(ic)*denw(ic) + &
                (1._kdp-phi(ic))*df(ic)*phfwt(ic)*dth(ic)
           f(ic) = pswp*phi(ic)*denw(ic) + phi(ic)*satnw(ic)*pdwp
           g(ic) = pswh*phi(ic)*denw(ic) + phi(ic)*satnw(ic)*pdwh
        END IF
     END IF
     IF(errflag) THEN
        WRITE (fustdout,9015)  &
             '***** Pressure-Enthalpy table error for cell no. ',ic,' *****'
        WRITE (fuclog,9015)  &
             '***** Pressure-Enthalpy table error for cell no. ',ic,' *****'
9015    FORMAT(a,i6,a)
     END IF
  END DO
! ... End of loop calculating properties and derivatives at nodes
  CALL avgprpty     ! ... compute the average density and viscosity at cell faces
END SUBROUTINE prpty
