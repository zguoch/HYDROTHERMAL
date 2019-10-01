SUBROUTINE upstre
  ! ... Purpose:  To determine which node is "upstream" of the cell boundary
  ! ...     between two nodes for both steam and water flow in each coordinate direction.
  ! ...     Compute the upsteam weighting factor which is applied to the relative
  ! ...     permeability and the enthalpy in the flow and transport equations.
  ! ...        The upstream weighting factor ranges from 0 to 1. 
  ! ...     Values less than 0.5 are for "+"
  ! ...     direction flow and values greater than 0.5 are for "-" direction
  ! ...     flow in the coordinate system.  
  ! ...     Flow direction is determined by the potential difference
  ! ...     between two nodes using the average fluid density at the cell face.
  ! ...     If the potential of the node to the "-" side of a cell face
  ! ...     is greater than the potential to the "+" side,
  ! ...     flow will be in the "+" coordinate direction.  If the two potentials differ
  ! ...     by more than a ratio factor (POTDIF, input), the upstream
  ! ...     weighting factor is set to 0 or 1 depending on the direction of flow for full
  ! ...     upstream weighting.
  ! ...     When the ratio of the potentials of two ajacent nodes is less
  ! ...     than POTDIF, the upsteam factor is scaled using a log function.
  ! ...     If the potentials are equal, the upsteam factor is 0.5 for centered-in-space
  ! ...     weighting.
  ! ...        When adjacent cells have different water saturations, 
  ! ...     it is sometimes necessary to adjust the upstream weighting
  ! ...     factor when one phase becomes immobile in one of the cells. In
  ! ...     laterally adjacent cells, as the relative permeability of a phase
  ! ...     approaches 0, the upstream weighting factor is gradually scaled to
  ! ...     match??? the relative permeability of the other phase. The same technique is used
  ! ...     for vertically adjacent cells when both phases are flowing the
  ! ...     same direction.  When there is counter-flow (water down and steam
  ! ...     up), then the upstream weighting for water is gradually scaled to
  ! ...     1 and the weighting factor for steam is scaled to 0.
  ! ...   Note: It is important to note that small changes in the upstream
  ! ...     weighting can significantly impact the numerical solution. Because
  ! ...     the derivatives of the upstream weighting have not yet been included
  ! ...     in the NR interation, the upstream weighting factors are held constant
  ! ...     after the third NR iteration unless there has been a phase change.
  ! ...        Left and bottom boundary cell faces are given a weighting factor of
  ! ...     -99, and right and top boundary cell faces are given a weighting factor
  ! ...     of -99. 
  ! ... In the air-water region, use upstream weighting exclusively.
  !
  USE machine_constants, ONLY: kdp
  USE parameters              !!$*** for flint
  USE control
  USE fdeq, ONLY: potdif, usx, uwx, usy, uwy, usz, uwz, xds, xdw, yds, ydw, zds, zdw
  USE mesh
!!$  USE parameters
  USE variables
  IMPLICIT NONE
  INTEGER :: i, ic, icup, icxm, icym, iczm, ik, ikxm, ikzm,  &
       j, k, mxf, mxl, mxm, mxp, myf, myl, mym, myp, mzf, mzl, mzm, mzp
  REAL(KIND=kdp) :: potratio, scalefac, rpval, uworig, usorig
  REAL(KIND=kdp) :: potratios, potratiow
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.7 $//$Date: 2006/10/03 23:42:11 $'
  !     ------------------------------------------------------------------
  IF (lnri > 3) RETURN      ! ... Do not adjust upstream weighting after 3rd N-R iteration
  ! ... Set the threshold value for difference ratio in relative permeability used to 
  ! ...     gradually scale the upstream weighting; presently 1%.
  rpval = 0.01_kdp
  ! ...         X direction upstream determination       
  DO  j = 1,ny
     DO  k = 1,nz
        ! ... Set weighting variable for first and last cell faces in X row
        mxf = (k-1)*nxx + 1 + (j-1)*nxx*nz
        mxl = k*nxx + (j-1)*nxx*nz
        usx(mxf) = 1._kdp
        uwx(mxf) = 1._kdp
        usx(mxl) = 0._kdp
        uwx(mxl) = 0._kdp
        IF (nx == 1) CYCLE
        DO  i = 2,nx
           ik = (k-1)*nx + i
           ikxm = ik - 1
           ! ... skip if cell (or the one to the left of it) is not in the active region
           IF (np(ik,j) == 0 .OR. np(ikxm,j) == 0) CYCLE
           ic = (k-1)*nx + i + (j-1)*nxz
           icxm = ic - 1
           mxm = (k-1)*nxx + i + (j-1)*nxx*nz
           ! ...               ----- water -----
           ! ... Calculate the water upstream weighting factor for the cell face
           ! ...     using the average water density at the face for
           ! ...     determining the ratio of the potentials.
           ! ... Set weighting to 0.5 if the water saturation is 0 in both cells.
           IF(satnw(ic) == 0.0_kdp .AND. satnw(icxm) == 0.0_kdp) THEN
              uwx(mxm) = 0.5_kdp
           ELSEIF(ind(ic) == 5 .OR. ind(icxm) == 5) THEN     ! ... upstream for air/water zone
              IF(p(icxm) > p(ic)) THEN
                 uwx(mxm) = 0._kdp     ! ... Positive direction flow
              ELSEIF(p(icxm) <= p(ic)) THEN
                 uwx(mxm) = 1._kdp     ! ... Negative direction flow
              ENDIF
           ELSE
              potratio = (p(ic)+xdw(mxm)*grav*zz(k))/(p(icxm)+xdw(mxm)*grav*zz(k))
              uwx(mxm) = 0.5_kdp/LOG10(1.0_kdp+potdif)*LOG10(potratio) + 0.5_kdp
              uwx(mxm) = MAX(uwx(mxm),0._kdp)
              uwx(mxm) = MIN(uwx(mxm),1._kdp)
           END IF
           ! ...              ----- steam -----
           ! ... Calculate the steam upstream weighting factor for the cell face
           ! ...      using the average steam density at the face for
           ! ...      determining the ratio of the potentials.
           ! ... Set weighting to 0.5 if the water saturation is 0 in both cells or 
           ! ...      one of the cells is in the air/water zone
           IF ((satnw(ic) == 1._kdp .AND. satnw(icxm) == 1._kdp) .OR.  &
                (ind(ic) == 5 .OR. ind(icxm) == 5)) THEN
              usx(mxm) = 0.5_kdp
           ELSE
              potratio = (p(ic)+xds(mxm)*grav*zz(k))/(p(icxm)+xds(mxm)*grav*zz(k))
              usx(mxm) = 0.5_kdp/LOG10(1.0_kdp+potdif)*LOG10(potratio) + 0.5_kdp
              usx(mxm) = MAX(usx(mxm),0._kdp)
              usx(mxm) = MIN(usx(mxm),1._kdp)
           END IF
           ! ... Adjust upsteam weighting for low saturation phase if necessary
           ! ...     When the saturations of adjacent cells are not equal, check
           ! ...     for possible adjustment of upstream weighting.
           ! ... No adjustment of air/water zone cell face weights
           IF(satnw(ic) /= satnw(icxm) .AND. ind(ic) /= 5 .AND. ind(icxm) /= 5) THEN
              ! ... Save the original value of upstream weighting
              uworig= uwx(mxm)
              usorig= usx(mxm)
              ! ... Determine the upstream cell
              IF (MAX(uwx(mxm),usx(mxm)) > 0.5_kdp) THEN
                 icup = ic
              ELSE
                 icup = icxm
              END IF
              ! ... Gradually scale upstream water when relative permeability to water 
              ! ...     is less than threshold (rpval)
              IF (rw(icup) < rpval) THEN
                 scalefac = 1.0_kdp - ((rpval - rw(icup))/rpval)     ! ... range (0-1)
                 uwx(mxm) = scalefac*uworig + (1.0_kdp - scalefac)*usorig
              END IF
              ! ... Gradually scale upstream water when relative permeability to steam
              ! ...     is less than threshold (rpval)
              IF (rs(icup) < rpval) THEN
                 scalefac = 1.0_kdp - ((rpval - rs(icup))/rpval)     ! ... range (0-1)
                 usx(mxm) = scalefac*usorig + (1.0_kdp - scalefac)*uworig
              END IF
           END IF
        END DO
     END DO
  END DO
  ! ...        Y direction upstream determination       
  DO  k = 1,nz
     DO  i = 1,nx
        ! ... Set weighting variable for first and last cell faces in Y row
        myf = (k-1)*i + i
        myl = myf + nxyz
        usy(myf) = 1._kdp
        uwy(myf) = 1._kdp
        usy(myl) = 0._kdp
        uwy(myl) = 0._kdp
        IF (ny == 1) CYCLE
        DO  j = 2, ny
           ik = (k-1)*nx + i
           ! ... skip if cell (or the one to the front of it) is not in the active region
           IF (np(ik,j) == 0 .OR. np(ik,j-1) == 0) CYCLE
           ic = (k-1)*nx + i + (j-1)*nxz
           icym = ic - nxz
           mym = ic
           ! ...             ----- water -----
           ! ... Calculate the water upstream weighting factor for the cell face
           ! ...     using the average water density at the face for
           ! ...     determining the ratio of the potentials.
           ! ... Set weighting to 0.5 if the water saturation is 0 in both cells or 
           ! ...      one of the cells is in the air/water region
           IF (satnw(ic) == 0.0_kdp .AND. satnw(icym) == 0.0_kdp) THEN
              uwy(mym) = 0.5_kdp
           ELSEIF(ind(ic) == 5 .OR. ind(icym) == 5) THEN
              IF(p(icym) > p(ic)) THEN
                 uwy(mym) = 0._kdp     ! ... Positive direction flow
              ELSEIF(p(icym) <= p(ic)) THEN
                 uwy(mym) = 1._kdp     ! ... Negative direction flow
              ENDIF
           ELSE
              potratio = (p(ic)+ydw(mym)*grav*zz(k))/(p(icym)+ydw(mym)*grav*zz(k))
              uwy(mym) = 0.5_kdp/LOG10(1.0_kdp+potdif)*LOG10(potratio) + 0.5_kdp
              uwy(mym) = MAX(uwy(mym),0._kdp)
              uwy(mym) = MIN(uwy(mym),1._kdp)
           END IF
           ! ...              ----- steam -----
           ! ... Calculate the steam upstream weighting factor for the cell face
           ! ...      using the average steam density at the face for
           ! ...      determining the ratio of the potentials.
           ! ... Set weighting to 0.5 if the water saturation is 0 in both cells or 
           ! ...      one of the cells is in the air/water region
           IF ((satnw(ic) == 1.0_kdp .AND. satnw(icym) == 1.0_kdp) .OR.  &
                (ind(ic) == 5 .OR. ind(icym) == 5)) THEN
              usy(mym) = 0.5_kdp
           ELSE
              potratio = (p(ic)+yds(mym)*grav*zz(k))/(p(icym)+yds(mym)*grav*zz(k))
              usy(mym) = 0.5_kdp/LOG10(1.0_kdp+potdif)*LOG10(potratio) + 0.5_kdp
              usy(mym) = MAX(usy(mym),0._kdp)
              usy(mym) = MIN(usy(mym),1._kdp)
           END IF
           ! ... Adjust upsteam weighting for low saturation phase if necessary
           ! ...     When the saturation of adjacent cells is not equal, check
           ! ...     for possible adjustment of upstream weighting.
           ! ... No adjustment of air-water region cell face weights
           IF (satnw(ic) /= satnw(icym) .AND. ind(ic) /= 5 .AND. ind(icym) /= 5) THEN
              ! ... Save the original value of upstream weighting
              uworig= uwy(mym)
              usorig= usy(mym)
              ! ... Determine the upstream cell
              IF (MAX(uwy(mym),usy(mym)) > 0.5_kdp) THEN
                 icup = ic
              ELSE
                 icup = icym
              END IF
              IF (rw(icup) < rpval) THEN        ! ... Gradually scale upstream water
                 scalefac = 1.0_kdp - ((rpval - rw(icup))/rpval)
                 uwy(mym) = scalefac*uworig + (1.0_kdp - scalefac)*usorig
              END IF
              IF (rs(icup) < rpval) THEN        ! ... Gradually scale upstream steam
                 scalefac = 1.0_kdp - ((rpval - rs(icup))/rpval)
                 usy(mym) = scalefac*usorig + (1.0_kdp - scalefac)*uworig
              END IF
           END IF
        END DO
     END DO
  END DO
  ! ...        Z direction upstream determination       
  DO  j = 1,ny
     DO  i = 1,nx
        ! ... Set weighting variable for first and last cell faces in Z column
        mzf = i + (j-1)*nzz*nx
        mzl = nz*nx + i + (j-1)*nzz*nx
        usz(mzf) = 1._kdp
        uwz(mzf) = 1._kdp
        usz(mzl) = 0._kdp
        uwz(mzl) = 0._kdp
        IF (nz == 1) CYCLE
        DO  k = 2, nz
           ik = (k-1)*nx + i
           ikzm = ik - nx
           ! ... Skip if cell (or the one below it) is not in the active region
           IF (np(ik,j) == 0 .OR. np(ikzm,j) == 0) CYCLE
           ic = (k-1)*nx + i + (j-1)*nxz
           iczm = ic - nx
           mzm = (k-1)*nx + i + (j-1)*nx*nzz
           ! ...             ----- water -----
           ! ... Calculate the water upstream weighting factor for the cell face
           ! ...     using the average water density at the face for
           ! ...     determining the ratio of the potentials.
           ! ... Set weighting to 0.5 if the water saturation is 0 in both cells or 
           ! ...      one of the cells is in the air/water region
           IF(satnw(ic) == 0.0_kdp .AND. satnw(iczm) == 0.0_kdp) THEN
              uwz(mzm) = 0.5_kdp
           ELSEIF(ind(ic) == 5 .OR. ind(iczm) == 5) THEN
              IF((p(ic)+zdw(mzm)*grav*zz(k)) > (p(iczm)+zdw(mzm)*grav*zz(k-1))) THEN
                 uwz(mzm) = 1._kdp     ! ... Negative direction flow
              ELSEIF((p(ic)+zdw(mzm)*grav*zz(k)) <= (p(iczm)+zdw(mzm)*grav*zz(k-1))) THEN
                 uwz(mzm) = 0._kdp     ! ... Positive direction flow
              ENDIF
           ELSE
              potratiow = (p(ic)+zdw(mzm)*grav*zz(k))/(p(iczm)+zdw(mzm)*grav*zz(k-1))
              uwz(mzm) = .5_kdp/LOG10(1.0_kdp+potdif)*LOG10(potratiow) + .5_kdp
              uwz(mzm) = MAX(uwz(mzm),0._kdp)
              uwz(mzm) = MIN(uwz(mzm),1._kdp)
           END IF
           ! ...              ----- steam -----
           ! ... Calculate the steam upstream weighting factor for the cell face
           ! ...      using the average steam density at the face for
           ! ...      determining the ratio of the potentials.
           ! ... Set weighting to 0.5 if the water saturation is 0 in both cells or 
           ! ...      one of the cells is in the air/water region
           IF ((satnw(ic) == 1.0_kdp .AND. satnw(iczm) == 1.0_kdp).OR.  &
                (ind(ic) == 5 .OR. ind(iczm) == 5)) THEN
              usz(mzm) = 0.5_kdp
           ELSE
              potratios = (p(ic)+zds(mzm)*grav*zz(k))/(p(iczm)+zds(mzm)*grav*zz(k-1))
              usz(mzm) = .5_kdp/LOG10(1.0_kdp+potdif)*LOG10(potratios) + .5_kdp
              usz(mzm) = MAX(usz(mzm),0._kdp)
              usz(mzm) = MIN(usz(mzm),1._kdp)
           END IF
           ! ... Adjust upsteam weighting for low saturation phase if necessary
           ! ...     When the saturation of adjacent cells is not equal, check
           ! ...     for possible adjustment of upstream weighting.
           ! ... No adjustment of air-water region cell face weights
           IF (satnw(ic) /= satnw(iczm) .AND. ind(ic) /= 5 .AND. ind(iczm) /= 5) THEN
              ! ... Save the original value of upstream weighting
              uworig= uwz(mzm)
              usorig= usz(mzm)
              ! ... Case of water and steam both flowing in the same direction
              IF((uwz(mzm) >= 0.5_kdp .AND. usz(mzm) >= 0.5_kdp) .OR.  &
                   (uwz(mzm) <= 0.5_kdp .AND. usz(mzm) <= 0.5_kdp)) THEN
                 ! ... Determine the upstream cell
                 IF (MAX(uwz(mzm),usz(mzm)) > 0.5_kdp) THEN
                    icup = ic
                 ELSE
                    icup = iczm
                 END IF
                 IF (rw(icup) < rpval) THEN            ! ... Gradually scale upstream water
                    scalefac = 1.0_kdp - ((rpval - rw(icup))/rpval)
                    uwz(mzm) = scalefac*uworig + (1.0_kdp - scalefac)*usorig
                 END IF
                 IF (rs(icup) < rpval) THEN            ! ... Gradually scale upstream steam
                    scalefac = 1.0_kdp - ((rpval - rs(icup))/rpval)
                    usz(mzm) = scalefac*usorig + (1.0_kdp - scalefac)*uworig
                 END IF
              ELSE
                 ! ... Case of two phase flow in opposite directions;
                 ! ...     steam flows up and water flows down
                 ! ... Gradually scale upstream water when rel.perm water < rpval
                 ! ...     from 0.5 for Rw=0 to calculated value for Rw=rpval ???
                 ! ...     from (1.0-Usz(mzm)) for Rw=0 to calculated value for Rw=rpval ???
                 IF (rw(ic) < rpval) THEN
                    scalefac = 1.0_kdp - ((rpval - rw(ic))/rpval)
                    uwz(mzm) = scalefac*uworig + (1.0_kdp - scalefac)
                 END IF
                 ! ... Gradually scale upstream steam when rel.perm steam < rpval
                 ! ...     from 0.5 for Rs=0 to calculated value for Rs=rpval
                 ! ...     from (1.0-Uwz(mzm)) for Rs=0 to calculated value for Rs=rpval
                 IF (rs(iczm) < rpval) THEN
                    scalefac = 1.0_kdp - (( rpval - rs(iczm))/rpval)
                    usz(mzm) = scalefac*usorig + (1.0_kdp - scalefac)*0.0_kdp
                 END IF
              END IF
           END IF
        END DO
     END DO
  END DO
  ! ... Cells at boundary of active region
  ! ...    cell faces on positive coordinate side are assigned UWX=-99
  ! ...    cell faces on negative coordinate side are assigned UWX=-99
  DO  j = 1,ny
     DO  k = 1,nz
        DO  i = 1,nx
           ik = (k-1)*nx + i
           ! ... Identify boundary cells
           IF (np(ik,j) /= 0) CYCLE
           ! ...           X direction
           mxm = (k-1)*nxx + i + (j-1)*nxx*nz
           mxp = mxm + 1
           usx(mxp) = -99._kdp
           uwx(mxp) = -99._kdp
           usx(mxm) = -99._kdp
           uwx(mxm) = -99._kdp
           ! ...           Y direction
           IF (ny > 1) THEN
              mym = (k-1)*nx + i + (j-1)*nxz
              myp = mym + nxz
              usy(myp) = -99._kdp
              uwy(myp) = -99._kdp
              usy(mym) = -99._kdp
              uwy(mym) = -99._kdp
           END IF
           ! ...           Z direction
           mzm = (k-1)*nx + i + (j-1)*nx*nzz
           mzp = mzm + nx
           usz(mzp) = -99._kdp
           uwz(mzp) = -99._kdp
           usz(mzm) = -99._kdp
           uwz(mzm) = -99._kdp
        END DO
     END DO
  END DO
END SUBROUTINE upstre
