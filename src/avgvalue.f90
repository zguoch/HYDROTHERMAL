SUBROUTINE avgvalue(ica, icb, adw, avw, avdw, pavdwp, pavdwh,  &
     avdgw, pavdgwp, pavdgwh, ads, avs, avds,  &
     pavdsp, pavdsh, avdgs, pavdgsp, pavdgsh)
  ! ... Purpose:  To compute the average density and viscosity of
  ! ...      steam and water and their derivatives w.r.t. P and H
  ! ...      at the cell boundary between two nodes.
  ! ... Method:  Compute the average P & H for the region between two nodes.
  ! ...      Then determine the density and viscosity for that point.  If
  ! ...      that point lies in the two phase region then determine the
  ! ...      density and viscosity for the coexisting steam and water for
  ! ...      that average P & H. When the two points being averaged lie on
  ! ...      different sides of the two-phase region or in the supercritical
  ! ...      region, then special techniques are applied to determine the
  ! ...      average properties.  These techniques attempt to make smooth
  ! ...      transitions, so that there are no marked jumps in density or
  ! ...      viscosity, when a different technique is applied.
  ! ...        Note: derivatives include 0.5 term to account for the
  ! ...            fact that the derivative is taken w.r.t. the P or H
  ! ...            at either one of the nodes and NOT at the cell face
  USE machine_constants, ONLY: kdp
  USE f_units
  USE parameters          !!$*** for flint
  USE control
!!$  USE parameters
  USE variables
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ica, icb
  REAL(KIND=kdp), INTENT(OUT) :: adw, avw, avdw, pavdwp, pavdwh, avdgw, pavdgwp, pavdgwh, &
       ads, avs, avds, pavdsp, pavdsh, avdgs, pavdgsp, pavdgsh
  ! ...   ICA  - IC number of node
  ! ...   ICB  - IC number of adjacent node in the x,y, or z direction
  ! ...   ADW  - average density of water in the two nodes
  ! ...   AVW  - average viscosity of water in the two nodes
  ! ...   AVDW  - average density/viscosity of water in the two nodes
  ! ...   PAVDWP  - derivative of AVDW w.r.t. P
  ! ...   PAVDWH  - derivative of AVDW w.r.t. H
  ! ...   AVDGW  - average g*den^2/vis of water in the two nodes
  ! ...   PAVDGWP  - derivative of AVDGW w.r.t. P
  ! ...   PAVDGWH  - derivative of AVDGW w.r.t. H
  ! ...   ADS  - average density of steam in the two nodes
  ! ...   AVS  - average viscosity of steam in the two nodes
  ! ...   AVDS  - average density/viscosity of steam in the two nodes
  ! ...   PAVDSP  - derivative of AVDS w.r.t. P
  ! ...   PAVDSH  - derivative of AVDS w.r.t. H
  ! ...   AVDGS  - average g*den^2/vis of steam in the two nodes
  ! ...   PAVDGSP  - derivative of AVDGS w.r.t. P
  ! ...   PAVDGSH  - derivative of AVDGS w.r.t. H
  !
  LOGICAL :: errflag
  REAL(KIND=kdp) :: pavg, havg, hw, hs, pdwh, pdwp, pvwh, pvwp, pdsh,  &
       pdsp, pvsh, pvsp, adwtmp, avwtmp, pdwptmp,  &
       pvwptmp, adstmp, avstmp, pdsptmp, pvsptmp, htemp
  REAL(KIND=kdp) :: relwtr, relstm, wtfac
  REAL(KIND=kdp) :: dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9
  REAL(KIND=kdp), PARAMETER :: pcrit=220.55e6_kdp
!!$          , hcrit=20.86e9_kdp
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.10 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  ! ... Initialize variables
  adw = 0._kdp
  avw = 1.0e20_kdp
  avdw = 0._kdp
  pavdwp = 0._kdp
  pavdwh = 0._kdp
  avdgw = 0._kdp
  pavdgwp = 0._kdp
  pavdgwh = 0._kdp
  ads = 0._kdp
  avs = 1.0e20_kdp
  avds = 0._kdp
  pavdsp = 0._kdp
  pavdsh = 0._kdp
  avdgs = 0._kdp
  pavdgsp = 0._kdp
  pavdgsh = 0._kdp
  errflag = .FALSE.
  ! ... Calculate average P and H
  pavg = 0.5_kdp*(p(ica)+p(icb))
  havg = 0.5_kdp*(h(ica)+h(icb))
  IF (pavg >= pcrit) THEN       ! ... average P&H in supercritical pressure region
     CALL tblookup(pavg, havg, adw, pdwh, pdwp,dum1, dum2, dum3, avw, pvwh, pvwp, errflag)
     IF(errflag) THEN
        WRITE(fustdout,'(2a,2i6)') 'ERROR - Problem calculating average water and steam properties; ',  &
             'Node numbers: ',ica,icb
        WRITE(fuclog,'(2a,2i6)') 'ERROR - Problem calculating average water and steam properties; ',  &
             'Node numbers: ',ica,icb
        ierr(75) = .TRUE.
        RETURN
     END IF
     pdsh = pdwh
     pdsp = pdwp
     pvsh = pvwh
     pvsp = pvwp
     ! ...   assign water and steam the same density and viscosity
     ads = adw
     avs = avw
     ! ...          calculate den/vis and g*den^2/vis terms
     avdw = adw/avw
     avds = avdw
     avdgw = grav*adw*adw/avw
     avdgs = avdgw
     ! ...   calculate derivative of density/viscosity at cell boundary
     ! ...          water
     pavdwh = 0.5_kdp*(pdwh/avw-adw*pvwh/(avw*avw))
     pavdwp = 0.5_kdp*(pdwp/avw-adw*pvwp/(avw*avw))
     pavdgwh = 0.5_kdp*grav*adw/avw*(2._kdp*pdwh-adw*pvwh/avw)
     pavdgwp = 0.5_kdp*grav*adw/avw*(2._kdp*pdwp-adw*pvwp/avw)
     ! ...          steam
     pavdsh = pavdwh
     pavdsp = pavdwp
     pavdgsh = pavdgwh
     pavdgsp = pavdgwp
!!$     RETURN
  ELSE        ! ...    average P&H in subcritical pressure region
     ! ...     calculate enthalpies of sat'd water (HW) and steam (HS) at pressure
     CALL phsbndy1(pavg, hw, hs)
     IF (havg >= hw .AND. havg <= hs) THEN     ! ...    average P&H in 2-phase region
        CALL phsbndy3(pavg, hw, hs, dum1, dum2, adw, pdwp, avw, pvwp,  &
             ads, pdsp, avs, pvsp, dum3, dum4, dum5)
        pdwh = 0._kdp
        pdsh = 0._kdp
        pvwh = 0._kdp
        pvsh = 0._kdp
        IF (ind(ica) == 1 .AND. ind(icb) == 1) THEN
           ! ...    both of the nodes are in the compressed water region 
           ! ...           (possible due to the curvature of the water satn boundary,
           ! ...           or due to suppressing a phase change this iteration)
           ! ...           Use the density and viscosity of water and steam as
           ! ...           calculated above.  However, use the average P and an H just
           ! ...           into the water field to calculate the derivatives of the
           ! ...           properties of water
           htemp = hw - 0.00001e9_kdp
           CALL tblookup(pavg, htemp, dum1, pdwh, pdwp,dum2, dum3, dum4, dum5, pvwh, pvwp, errflag)
           IF(errflag) THEN
              WRITE(fustdout,'(2a,2i6)') 'ERROR - Problem calculating average water and steam properties; ',  &
                   'Node numbers: ',ica,icb
              WRITE(fuclog,'(2a,2i6)') 'ERROR - Problem calculating average water and steam properties; ',  &
                   'Node numbers: ',ica,icb
              ierr(75) = .TRUE.
              RETURN
           END IF
        ELSE IF (ind(ica) == 3 .AND. ind(icb) == 3) THEN
           ! ...    both of the nodes are in the superheated steam region 
           ! ...           (possible due to the curvature of the steam satn boundary,
           ! ...           or due to suppressing a phase change this iteration)
           ! ...           Use the density and viscosity of water and steam as
           ! ...           calculated above.  However, use the average P and an H just
           ! ...           into the steam region to calculate the derivatives of the
           ! ...           properties of steam
           htemp = hs + 0.00001e9_kdp
           CALL tblookup(pavg, htemp, ads, pdsh, pdsp,dum1, dum2, dum3, avs, pvsh, pvsp, errflag)
           IF(errflag) THEN
              WRITE(fustdout,'(2a,2i6)') 'ERROR - Problem calculating average water and steam properties; ',  &
                   'Node numbers: ',ica,icb
              WRITE(fuclog,'(2a,2i6)') 'ERROR - Problem calculating average water and steam properties; ',  &
                   'Node numbers: ',ica,icb
              ierr(75) = .TRUE.
              RETURN
           END IF
           pdwp = 0._kdp
           pdwh = 0._kdp
           pvwp = 0._kdp
           pvwh = 0._kdp
           adw = 0._kdp
           avw = 1.0e20_kdp
        END IF
        ! ...          calculate den/vis and g*den^2/vis terms
        avdw = adw/avw
        avds = ads/avs
        avdgw = grav*adw*adw/avw
        avdgs = grav*ads*ads/avs
        ! ...        calculate derivative of density/viscosity at cell boundary
        ! ...             water
        pavdwh = 0._kdp
        pavdwp = 0.5_kdp*(pdwp/avw-adw*pvwp/(avw*avw))
        pavdgwh = 0._kdp
        pavdgwp = 0.5_kdp*grav*adw/avw*(2._kdp*pdwp-adw*pvwp/avw)
        ! ...             steam
        pavdsh = 0._kdp
        pavdsp = 0.5_kdp*(pdsp/avs-ads*pvsp/(avs*avs))
        pavdgsh = 0._kdp
        pavdgsp = 0.5_kdp*grav*ads/avs*(2._kdp*pdsp-ads*pvsp/avs)
!!$        RETURN
     ELSE IF (havg < hw) THEN        ! ... average P&H in compressed water region
        CALL tblookup(pavg, havg, adw, pdwh, pdwp,dum1, dum2, dum3, avw, pvwh, pvwp, errflag)
        IF(errflag) THEN
           WRITE(fustdout,'(2a,2i6)') 'ERROR - Problem calculating average water and steam properties; ',  &
                'Node numbers: ',ica,icb
           WRITE(fuclog,'(2a,2i6)') 'ERROR - Problem calculating average water and steam properties; ',  &
                'Node numbers: ',ica,icb
           ierr(75) = .TRUE.
           RETURN
        END IF
        IF ((ind(ica) == 1 .OR. ind(ica) == 5) .AND. p(ica) <= pcrit .AND.  &
             (ind(icb) == 1 .OR. ind(icb) == 5) .AND. p(icb) <= pcrit) THEN
           ! ...   both nodes are in compressed water region or in the air-water region
           ! ...            Make water and steam values identical if one node has
           ! ...            P>220.0   (as is done in prpty)
           IF (p(ica) > 220.0e6_kdp .OR. p(icb) > 220.0e6_kdp) THEN
              ads = adw
              avs = avw
              pdsp = pdwp
              pdsh = pdwh
              pvsp = pvwp
              pvsh = pvwh
           ELSE
              ads = 0._kdp
              avs = 1.0e20_kdp
              pdsp = 0._kdp
              pdsh = 0._kdp
              pvsp = 0._kdp
              pvsh = 0._kdp
           END IF
        ELSE IF (ind(ica) == 2 .AND. ind(icb) == 2) THEN
           ! ...    both of the nodes are in the 2-phase region 
           ! ...           (possible due to the curvature of the water satn boundary,
           ! ...           or due to suppressing a phase change this iteration).
           ! ...           Use the average P to calculate properties of steam.
           ! ...           For water, use the density and viscosity determined from
           ! ...           the avg P&H, but use the derivatives for water from
           ! ...           along satn curve
           CALL phsbndy3(pavg, hw, hs, dum1, dum2, dum3, pdwp, dum4, pvwp,  &
                ads, pdsp, avs, pvsp, dum5, dum6, dum7)
           pdwh = 0._kdp
           pvwh = 0._kdp
           pdsh = 0._kdp
           pvsh = 0._kdp
        ELSE IF (ind(ica) == 2 .OR. ind(icb) == 2 .OR.  &
             ind(ica) == 3 .OR. ind(icb) == 3) THEN
           ! ...    one of the nodes is in the 2-phase region or in the
           ! ...           superheated steam region.
           ! ...           Calculate the density and
           ! ...           viscosity of a hypothetical steam using the average pressure
           CALL phsbndy3(pavg, hw, hs, dum1, dum2, dum3, dum4, dum5, dum6,  &
                ads, pdsp, avs, pvsp, dum7, dum8, dum9)
           pdsh = 0._kdp
           pvsh = 0._kdp
        ELSE IF (ind(ica) == 4 .OR. ind(icb) == 4) THEN
           ! ...    one of the nodes is in the supercritical high H region
           ! ...            Calculate the density and viscosity of a hypothetical
           ! ...            steam using the average pressure, and then use a weighted
           ! ...            average between the water properties and the hypothetical
           ! ...            steam, based on the gradiational relative permeability
           ! ...            (determined in relperm2) for the endpoint in the super
           ! ...            critical region
           CALL phsbndy3(pavg, hw, hs, dum1, dum2, dum3, dum4, dum5, dum6,  &
                adstmp, pdsptmp, avstmp, pvsptmp, dum7, dum8, dum9)
           IF (ind(ica) == 4) THEN
              CALL relperm2(p(ica), h(ica), dum1, dum2, dum3, relstm, dum4, dum5)
           ELSE
              CALL relperm2(p(icb), h(icb), dum1, dum2, dum3, relstm, dum4, dum5)
           END IF
           ! ...           relative perm of steam in this region should be 0.5 - 1.0
           wtfac = (relstm-0.5_kdp)*2._kdp
           pdsp = wtfac*pdsptmp + (1._kdp-wtfac)*pdwp
           pdsh = wtfac*0._kdp + (1._kdp-wtfac)*pdwh
           pvsp = wtfac*pvsptmp + (1._kdp-wtfac)*pvwp
           pvsh = wtfac*0._kdp + (1._kdp-wtfac)*pvwh
           ads = wtfac*adstmp + (1._kdp-wtfac)*adw
           avs = wtfac*avstmp + (1._kdp-wtfac)*avw
        ELSE IF ((ind(ica) == 1 .AND. p(ica) > pcrit) .OR.  &
             (ind(icb) == 1 .AND. p(icb) > pcrit)) THEN
           ! ...    one of the nodes is in the supercritical low H region
           ! ...           Assign water and steam the same density and viscosity
           pdsh = pdwh
           pdsp = pdwp
           pvsh = pvwh
           pvsp = pvwp
           ads = adw
           avs = avw
        ELSE
           WRITE(fustdout,*) 'ERROR - Problem calculating average steam properties'
           ierr(121) = .TRUE.
           RETURN
        END IF
        avdw = adw/avw
        avds = ads/avs
        avdgw = grav*adw*adw/avw
        avdgs = grav*ads*ads/avs
        ! ...        calculate derivative of density/viscosity at cell boundary
        ! ...             water
        pavdwh = 0.5_kdp*(pdwh/avw-adw*pvwh/(avw*avw))
        pavdwp = 0.5_kdp*(pdwp/avw-adw*pvwp/(avw*avw))
        pavdgwh = 0.5_kdp*grav*adw/avw*(2._kdp*pdwh-adw*pvwh/avw)
        pavdgwp = 0.5_kdp*grav*adw/avw*(2._kdp*pdwp-adw*pvwp/avw)
        ! ...             steam
        pavdsh = 0.5_kdp*(pdsh/avs-ads*pvsh/(avs*avs))
        pavdsp = 0.5_kdp*(pdsp/avs-ads*pvsp/(avs*avs))
        pavdgsh = 0.5_kdp*grav*ads/avs*(2._kdp*pdsh-ads*pvsh/avs)
        pavdgsp = 0.5_kdp*grav*ads/avs*(2._kdp*pdsp-ads*pvsp/avs)
!!$        RETURN
     ELSE IF (havg > hs) THEN        ! ...  for average P&H in superheated steam region
        CALL tblookup(pavg, havg, ads, pdsh, pdsp,dum1, dum2, dum3, avs, pvsh, pvsp, errflag)
        IF(errflag) THEN
           WRITE(fustdout,'(2a,2i6)') 'ERROR - Problem calculating average water and steam properties; ',  &
                'Node numbers: ',ica,icb
           WRITE(fuclog,'(2a,2i6)') 'ERROR - Problem calculating average water and steam properties; ',  &
                'Node numbers: ',ica,icb
           ierr(75) = .TRUE.
           RETURN
        END IF
        IF(ind(ica) == 3 .AND. ind(icb) == 3) THEN
           ! ...  both nodes are in superheated steam region. Make water and steam values 
           ! ...      identical if one node has P>220.e6   (as is done in prpty)
           IF (p(ica) > 220.e6_kdp .OR. p(icb) > 220.e6_kdp) THEN
              adw = ads
              avw = avs
              pdwp = pdsp
              pdwh = pdsh
              pvwp = pvsp
              pvwh = pvsh
           ELSE
              adw = 0._kdp
              avw = 1.0e20_kdp
              pdwp = 0._kdp
              pdwh = 0._kdp
              pvwp = 0._kdp
              pvwh = 0._kdp
           END IF
        ELSE IF (ind(ica) == 2 .AND. ind(icb) == 2) THEN
           ! ...    both of the nodes are in the 2-phase region 
           ! ...           (possible due to the curvature of the steam satn boundary,
           ! ...           or due to suppressing a phase change this iteration)
           ! ...           Use the average P to calculate properties of water.
           ! ...           For steam, use the density and viscosity determined from
           ! ...           the avg P&H, but use the derivatives for steam from
           ! ...           along satn curve
           CALL phsbndy3(pavg, hw, hs, dum1, dum2, adw, pdwp, avw, pvwp,  &
                dum3, pdsp, dum4, pvsp, dum5, dum6, dum7)
           pdwh = 0._kdp
           pvwh = 0._kdp
           pdsh = 0._kdp
           pvsh = 0._kdp
        ELSE IF (ind(ica) == 2 .OR. ind(icb) == 2 .OR.  &
             (ind(ica) == 1 .AND. p(ica) <= pcrit) .OR.  &
             (ind(icb) == 1 .AND. p(icb) <= pcrit)) THEN
           ! ...   one of the nodes is in the 2-phase region or in the
           ! ...           compressed water region. Calculate the density and
           ! ...           viscosity of a hypothetical water using the average pressure
           CALL phsbndy3(pavg, hw, hs, dum1, dum2, adw, pdwp, avw, pvwp,  &
                dum3, dum4, dum5, dum6, dum7, dum8, dum9)
           pdwh = 0._kdp
           pvwh = 0._kdp
        ELSE IF ((ind(ica) == 1 .AND. p(ica) > pcrit) .OR.  &
             (ind(icb) == 1 .AND. p(icb) > pcrit)) THEN
           ! ...    one of the nodes is in the supercritical low H region
           ! ...            Calculate the density and viscosity of a hypothetical
           ! ...            water using the average pressure, and then use a weighted
           ! ...            average between the steam properties and the hypothetical
           ! ...            water, based on the gradiational relative permeability
           ! ...            (determined in relperm2) for the endpoint in the super
           ! ...            critical region
           CALL phsbndy3(pavg, hw, hs, dum1, dum2,  &
                adwtmp, pdwptmp, avwtmp, pvwptmp, dum3, dum4, dum5, dum6,  &
                dum7, dum8, dum9)
           IF (ind(ica) == 1 .AND. p(ica) > pcrit) THEN
              CALL relperm2(p(ica), h(ica), relwtr, dum1, dum2, dum3, dum4, dum5)
              CALL relperm2(p(icb), h(icb), relwtr, dum1, dum2, dum3, dum4, dum5)
           END IF
           ! ...           relative perm of water in this region should be 0.0 - 0.5
           wtfac = relwtr*2._kdp
           pdwp = wtfac*pdwptmp + (1._kdp-wtfac)*pdsp
           pdwh = wtfac*0._kdp + (1._kdp-wtfac)*pdsh
           pvwp = wtfac*pvwptmp + (1._kdp-wtfac)*pvsp
           pvwh = wtfac*0._kdp + (1._kdp-wtfac)*pvsh
           adw = wtfac*adwtmp + (1._kdp-wtfac)*ads
           avw = wtfac*avwtmp + (1._kdp-wtfac)*avs
        ELSE IF (ind(ica) == 4 .OR. ind(icb) == 4) THEN
           ! ...   one of the nodes is in the supercritical high H region
           ! ...           Assign water and steam the same density and viscosity
           pdwh = pdsh
           pdwp = pdsp
           pvwh = pvsh
           pvwp = pvsp
           adw = ads
           avw = avs
        ELSE
           WRITE (fustdout,*) 'ERROR - Problem calculating average water properties'
           ierr(74) = .TRUE.
           RETURN
        END IF
        avdw = adw/avw
        avds = ads/avs
        avdgw = grav*adw*adw/avw
        avdgs = grav*ads*ads/avs
        ! ...        calculate derivative of density/viscosity at cell boundary
        ! ...             water
        pavdwh = 0.5_kdp*(pdwh/avw-adw*pvwh/(avw*avw))
        pavdwp = 0.5_kdp*(pdwp/avw-adw*pvwp/(avw*avw))
        pavdgwh = 0.5_kdp*grav*adw/avw*(2._kdp*pdwh-adw*pvwh/avw)
        pavdgwp = 0.5_kdp*grav*adw/avw*(2._kdp*pdwp-adw*pvwp/avw)
        ! ...             steam
        pavdsh = 0.5_kdp*(pdsh/avs-ads*pvsh/(avs*avs))
        pavdsp = 0.5_kdp*(pdsp/avs-ads*pvsp/(avs*avs))
        pavdgsh = 0.5_kdp*grav*ads/avs*(2._kdp*pdsh-ads*pvsh/avs)
        pavdgsp = 0.5_kdp*grav*ads/avs*(2._kdp*pdsp-ads*pvsp/avs)
!!$        RETURN
     END IF
  END IF
END SUBROUTINE avgvalue
