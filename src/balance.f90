SUBROUTINE balance
  ! ... Purpose:  Computes the mass and energy balance tables for the simulation.
  ! ...        The mass and energy contributions from the constant value nodes and
  ! ...        seepage nodes and
  ! ...        the source nodes are determined.  The change in mass and
  ! ...        energy in storage is also calculated.  
  ! ...        The residual mass and energy represents the amount of imbalance in the
  ! ...        flow or energy transport finite difference equations.
  ! ...        This is reported along with the percent change in mass and energy
  ! ...        in the region.
  USE machine_constants, ONLY: kdp
  USE f_units
  USE bc
  USE control
  USE mesh
  USE source
  USE units
  USE variables
  IMPLICIT NONE
  INTEGER :: i, ic, ik, j, k
  LOGICAL :: lprint
  REAL(KIND=kdp) :: ahesumi, ahesumo, cheri, chero, chesumi, chesumo, cpmri, cpmro,  &
       cpmsumi, cpmsumo, dehir,  &
       dfir, eeei, eeeo, eeeri, eeero, ehir, ensum, gveqne, gveqnm, fir, hfesumi, &
       lqir,  & 
       pfmsumi, pfesumi,  &
       qqqi, qqqo,  &
       qqqri, qqqro,  &
       qsppbc, qseepbc, &
       scrfir,  &
       seepmri, seepmro, seeperi, seepero,  &
       seepmsumi, seepmsumo, seepesumi, seepesumo,  &
       sfres, sfrespct, shres, shrespct, stotfi, stotfp, stothi, stothp,  &
       tdfir, tdehir, tfrespct, threspct,  &
       xmsum, u1, upri, vpir
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.16 $//$Date: 2007/10/16 21:17:22 $'
  !     ------------------------------------------------------------------
  !...
  stotfi = 0._kdp
  stotfp = 0._kdp
  stothi = 0._kdp
  stothp = 0._kdp
  ! ... Compute the total mass and energy contributions from the
  ! ...      specified pressure and enthalpy b.c. nodes
  ! ...      and the specified pressure with associated enthalpy b.c. nodes
  cpmsumi = 0._kdp
  cpmsumo = 0._kdp
  chesumi = 0._kdp
  chesumo = 0._kdp
  ahesumi = 0._kdp
  ahesumo = 0._kdp
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF (ibc(ik,j) == 11) THEN          ! ... Select specified P and H nodes
        CALL goveqn(gveqnm, gveqne, i, j, k)
        ! ... The total mass and energy contribution for the current time step is
        ! ... the negative of the governing equation times the time step length
        ! ... because change in storage is zero and no sources are allowed in a specified
        ! ...      P and H cell
        IF(gveqnm < 0._kdp) THEN
           cpmsumi = cpmsumi - gveqnm*delt
        ELSE
           cpmsumo = cpmsumo + gveqnm*delt
        END IF
        IF(gveqne < 0._kdp) THEN
           chesumi = chesumi - gveqne*delt
        ELSE
           chesumo = chesumo + gveqne*delt
        END IF
     ELSEIF (ibc(ik,j) == 10) THEN          ! ... Select specified P and associated H nodes
        CALL gov_spp(gveqnm,i,j,k)
        ! ... Add in the capacitance and source terms
        gveqnm = gveqnm - (xm(ic)-xmoldt(ic))/delt + q(ic)
        qsppbc = -gveqnm      ! ... Change sign so inflow is positive
        IF(qsppbc > 0._kdp) THEN     ! ... inflow
           gveqne = qsppbc*ehassoc(ic)
           cpmsumi = cpmsumi + qsppbc*delt
           ahesumi = ahesumi + gveqne*delt
        ELSE                         ! ... outflow
           gveqne = qsppbc*h(ic)
           cpmsumo = cpmsumo - qsppbc*delt
           ahesumo = ahesumo - gveqne*delt
        END IF
     END IF
  END DO
  tcfsbc = tcfsbc + cpmsumi - cpmsumo
  tchsbc = tchsbc + chesumi - chesumo + ahesumi - ahesumo
  stotfi = stotfi + cpmsumi
  stotfp = stotfp + cpmsumo
  stothi = stothi + chesumi + ahesumi
  stothp = stothp + chesumo + ahesumo
  totfi = totfi + cpmsumi
  totfp = totfp + cpmsumo
  tothi = tothi + chesumi + ahesumi
  tothp = tothp + chesumo + ahesumo
  ! ... Compute the total mass and energy contributions from the seepage cells
  seepmsumi = 0._kdp
  seepmsumo = 0._kdp
  seepesumi = 0._kdp
  seepesumo = 0._kdp
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     ! ... Select seepage cells that are seeping
     IF ((ibc(ik,j) == 30 .OR. ibc(ik,j) == 50) .AND. seeping(ik,j)) THEN
        CALL gov_seep(gveqnm,gveqne,i,j,k)
        ! ... Add in the capacitance term for flow.
        ! ... No precipitation source term for flow. Problem if a well is
        ! ...     also in a seeping cell.
        gveqnm = gveqnm - (xm(ic)-xmoldt(ic))/delt
        qseepbc = -gveqnm      ! ... Change sign so inflow is positive        
        ! ... The seepage mass contribution for the current time step is
        ! ... the negative of the governing equation times the time step length
        ! ... The seepage energy contribution for the current time step is
        ! ... the advective heat outflow from the flow equation.
        IF(qseepbc < 0._kdp) THEN    ! ... outflow
           gveqne = qseepbc*h(ic)
           seepmsumo = seepmsumo - qseepbc*delt
           seepesumo = seepesumo - gveqne*delt
        END IF
     END IF
  END DO
  tcfsepbc = tcfsepbc + seepmsumi - seepmsumo
  tchsepbc = tchsepbc + seepesumi - seepesumo
  stotfi = stotfi + seepmsumi
  stotfp = stotfp + seepmsumo
  stothi = stothi + seepesumi
  stothp = stothp + seepesumo
  totfi = totfi + seepmsumi
  totfp = totfp + seepmsumo
  tothi = tothi + seepesumi
  tothp = tothp + seepesumo
  ! ... Compute the total mass and energy contributions from the precipitation flux cells
  pfmsumi = 0._kdp
  pfesumi = 0._kdp
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF (ibc(ik,j) == 20) THEN          ! ... Select precipitation flux cells
        pfmsumi = pfmsumi + q(ic)*delt
        pfesumi = pfesumi + qh(ic)*delt
     ELSEIF (ibc(ik,j) == 50 .AND. .NOT.seeping(ik,j)) THEN   
        ! ... Select precipitation recharge flux cells where seepage is off
        pfmsumi = pfmsumi + q(ic)*delt
        pfesumi = pfesumi + qh(ic)*delt
     END IF
  END DO
  tcfpfbc = tcfpfbc + pfmsumi
  tchpfbc = tchpfbc + pfesumi
  stotfi = stotfi + pfmsumi
  stothi = stothi + pfesumi
  totfi = totfi + pfmsumi
  tothi = tothi + pfesumi
  ! ... Compute the total energy contributions from the basal heat flux cells
  hfesumi = 0._kdp
  DO ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF (ibc(ik,j) == 02) THEN          ! ... Select basal heat flux cells
        hfesumi = hfesumi + qhcondflx(ic)*delt
     END IF
  END DO
  tchhcbc = tchhcbc + hfesumi
  stothi = stothi + hfesumi
  tothi = tothi + hfesumi
  ! ... Compute the mass and energy contributions from the source nodes;
  ! ...      node points and well layers
  qqqi = 0._kdp
  qqqo = 0._kdp
  eeei = 0._kdp
  eeeo = 0._kdp
  DO ic=1,nxyz
     IF(iq(ic) == 1) THEN                ! ... Select source cells
        IF(qnodsrc(ic) > 0._kdp) THEN
           qqqi = qqqi + qnodsrc(ic)*delt
        ELSE
           qqqo = qqqo - qnodsrc(ic)*delt
        END IF
        IF(qhnodsrc(ic) > 0._kdp) THEN
           eeei = eeei + qhnodsrc(ic)*delt
        ELSE
           eeeo = eeeo - qhnodsrc(ic)*delt
        END IF
        IF(qwelsrc(ic) > 0._kdp) THEN
           qqqi = qqqi + qwelsrc(ic)*delt
        ELSE
           qqqo = qqqo - qwelsrc(ic)*delt
        END IF
        IF(qhwelsrc(ic) > 0._kdp) THEN
           eeei = eeei + qhwelsrc(ic)*delt
        ELSE
           eeeo = eeeo - qhwelsrc(ic)*delt
        END IF
     END IF
  END DO
  totwfi = totwfi + qqqi
  totwfp = totwfp + qqqo
  totwhi = totwhi + eeei
  totwhp = totwhp + eeeo
  stotfi = stotfi + qqqi
  stotfp = stotfp + qqqo
  stothi = stothi + eeei
  stothp = stothp + eeeo
  totfi = totfi + qqqi
  totfp = totfp + qqqo
  tothi = tothi + eeei
  tothp = tothp + eeeo
  ! ... Compute the change in mass and energy storage and
  ! ...    sum the total mass and energy in the simulation region
  xmsum = 0._kdp
  ensum = 0._kdp
  dfir = 0._kdp
  dehir = 0._kdp
  DO  ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     IF (ibc(ik,j) /= -1) THEN        ! ... Select active cells
        xmsum = xmsum + xm(ic)
        ensum = ensum + en(ic)
        dfir = dfir + (xm(ic)-xmoldt(ic))
        dehir = dehir + (en(ic)-enoldt(ic))
     END IF
  END DO
  ! ... Calculate amount of vapor, liquid, and super-critical fluid in system
  CALL steam_mass(vpir, lqir, scrfir)
  fir = xmsum
  ehir = ensum
  sfres = dfir - stotfi + stotfp
  tfres = fir - fir0 - totfi + totfp
  u1 = MAX(ABS(dfir),stotfi,stotfp)
  IF(u1 > 0.) sfrespct = 100._kdp*sfres/u1
  u1 = MAX(ABS(fir-fir0),totfi,totfp)
  IF(u1 > 0.) tfrespct = 100._kdp*tfres/u1
  shres = dehir - stothi + stothp
  thres = ehir - ehir0 - tothi + tothp
  shrespct = 100._kdp*shres/MAX(ABS(dehir),stothi,stothp)
  threspct = 100._kdp*thres/MAX(ABS(ehir-ehir0),tothi,tothp)
  ! ... Compute the rates of change
  ! **** not used for printout at present ****
  ! ... Specified pressure and enthalpy
  cpmri = cpmsumi/delt
  cpmro = cpmsumo/delt
  cheri = chesumi/delt
  chero = chesumo/delt
  ! ... Seepage
  seepmri = seepmsumi/delt
  seepmro = seepmsumo/delt
  seeperi = seepesumi/delt
  seepero = seepesumo/delt
  ! ... Heat flux

  ! ... Precipitation flux

  ! ... Source/sink node points and wells
  qqqri = qqqi/delt
  qqqro = qqqo/delt
  eeeri = eeei/delt
  eeero = eeeo/delt
  IF(ABS(print_balance) > 0.) THEN
     CALL print_control(print_balance,time,ltime,tchg,timprbal,lprint)
     IF (lprint) THEN
        WRITE(fubal,2001)  '*** Output at End of Time Step No. ',ltime,' ***'
2001    FORMAT(/tr20,a,i5,a)
        IF (iyr == 2) THEN
           WRITE(fubal,2002) 'Time '//dots,year(time),'(yr)'
2002       FORMAT(/tr5,a52,1pg12.3,tr1,a)
        ELSE
           WRITE(fubal,2002) 'Time '//dots,time,'(s)'
        END IF
        WRITE(fubal,2010)
2010    FORMAT(/tr30,'*** Global Flow Balance Summary ***'/tr15,  &
             'Current Time Step',tr25,'Rates',tr18,'Amounts')
        WRITE(fubal,2011) 'Fluid inflow '//dots,stotfi/delt,'(g/s)',stotfi,'(g)',  &
             'Fluid outflow '//dots,stotfp/delt,'(g/s)',stotfp,'(g)',  &
             'Change in fluid in region '//dots,dfir/delt,'(g/s)',dfir,'(g)',  &
             'Residual imbalance '//dots,sfres/delt,'(g/s)',sfres,'(g)',  &
             'Percent imbalance '//dots,sfrespct
2011    FORMAT(/4(tr1,a52,1pe14.6,tr2,a,tr3,e14.6,tr1,a/),tr1,a52,tr28,0pf8.2)
        WRITE(fubal,2011) 'Heat inflow '//dots,stothi/delt,'(erg/s)',stothi,'(erg)',  &
             'Heat outflow '//dots,stothp/delt,'(erg/s)',stothp,'(erg)',  &
             'Change in heat in region '//dots,dehir/delt,'(erg/s)',dehir,'(erg)',  &
             'Residual imbalance '//dots,shres/delt,'(erg/s)',shres,'(erg)',  &
             'Percent imbalance '//dots,shrespct
        WRITE(fubal,2015) 'Amounts'
2015    FORMAT(/tr55,a)
        WRITE(fubal,2023) 'Step total specified p cell fluid net inflow '//  &
             dots,cpmsumi-cpmsumo,'(g)',  &
             'Step total precipitation flux b.c. fluid net inflow '//dots,pfmsumi,'(g)',  &
             'Step total seepage b.c. fluid net inflow '//dots,seepmsumi-seepmsumo,'(g)',  &
             'Step total source fluid net inflow '//dots,qqqi-qqqo, '(g)'
2023    FORMAT(/6(tr1,a52,1pe14.6,tr1,a/))
        WRITE(fubal,2123) 'Step total specified h (or T) or associated with',  &
             '   specified p cell heat net inflow '//dots,  &
             chesumi-chesumo+ahesumi-ahesumo,'(erg)',  &
             'Step total precipitation flux b.c. heat net inflow '//dots,pfesumi,'(erg)',  &
             'Step total seepage b.c. heat net inflow '//dots,seepesumi-seepesumo,'(erg)',  &
             'Step total heat flux b.c. heat net inflow '//dots,hfesumi,'(erg)',  &
             'Step total source heat net inflow '//dots,eeei-eeeo,'(erg)'
2123    FORMAT(/tr1,a/6(tr1,a52,1pe14.6,tr1,a/))
        WRITE(fubal,2017) 'Cumulative Summary','Amounts'
2017    FORMAT(/tr15,a/tr55,a)
        tdfir=fir-fir0
        WRITE(fubal,2018) 'Fluid inflow '//dots,totfi,'(g)',  &
             'Fluid outflow '//dots,totfp,'(g)',  &
             'Change in fluid in region '//dots,tdfir,'(g)',  &
             'Fluid (liquid + vapor + super-critical) in region '//dots,fir,'(g)',  &
             '   Liquid in region '//dots,lqir,'(g)',  &
             '   Vapor in region '//dots,vpir,'(g)',  &
             '   Super-critical fluid in region '//dots,scrfir,'(g)',  &
             'Residual imbalance '//dots,tfres,'(g)',  &
             'Percent imbalance '//dots,tfrespct
2018    FORMAT(/8(tr1,a52,1pe14.6,tr1,a/),tr1,a52,0pf12.2)
        tdehir=ehir-ehir0
        WRITE(fubal,2019) 'Heat inflow '//dots,tothi,'(erg)',  &
             'Heat outflow '//dots,tothp,'(erg)',  &
             'Change in heat in region '//dots,tdehir,'(erg)',  &
             'Heat in region '//dots,ehir,'(erg)',  &
             'Residual imbalance '//dots,thres,'(erg)',  &
             'Percent imbalance '//dots,threspct
2019    FORMAT(/5(tr1,a52,1pe14.6,tr1,a/),tr1,a52,0pf12.2)
        WRITE(fubal,2023) 'Cumulative specified p cell fluid net inflow '//  &
             dots,tcfsbc,'(g)',  &
             'Cumulative precipitation flux b.c. fluid net inflow '//dots,tcfpfbc,'(g)',  &
             'Cumulative seepage b.c. fluid net inflow '//dots,tcfsepbc,'(g)',  &
             'Cumulative source fluid net inflow '//dots,totwfi-totwfp,'(g)'
        WRITE(fubal,2123) 'Cumulative specified h (or T) or associated with',  &
             '   specified p cell heat net inflow '//  &
             dots,tchsbc,'(erg)',  &
             'Cumulative precipitation flux b.c. heat net inflow '//dots,tchpfbc,'(erg)',  &
             'Cumulative seepage b.c. heat net inflow '//dots,tchsepbc,'(erg)',  &
             'Cumulative heat flux b.c. heat net inflow '//dots,tchhcbc,'(erg)',  &
             'Cumulative source fluid heat net inflow '//dots,totwhi-totwhp,'(erg)'
     END IF
     ! ... Set the next time for printout
     prtimenxt = MIN(tchg, timprp, timprh, timprt, timprsat, timprd, timprv, timprpot,  &
          timprvel, timprbcf, timprsrc, timprprop, timprpor, timprperm, timprbal, timprdump,  &
          timprresid, timprdimno, timprpscal, timprpvec, timprtemp)
  END IF
END SUBROUTINE balance
